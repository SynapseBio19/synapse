"""ORF search benchmarked against the RefSeq gene annotation of E. coli K-12.

Reference data
    *Escherichia coli* K-12 substr. MG1655, RefSeq NC_000913.3: genome
    (FASTA) and its annotation (GFF3) from NCBI, cached in
    ``.cache/benchmarks/`` and verified by SHA-256. The annotated CDS are
    the reference ("truth set"); they are curated gene models, not ORFs.

Method
    Synapse searches all six frames with the bacterial code (table 11),
    the start codons ATG, GTG and TTG, no nesting (one ORF per stop: the
    longest, i.e. the most upstream start) and a 90-nt minimum, the
    smallest length that every evaluated CDS reaches except a few leader
    peptides (they are counted as misses, not dropped).

    A CDS is evaluated when it is a single-segment, non-pseudo CDS
    without a ribosomal-slippage exception. Selenoprotein genes (internal
    UGA read as selenocysteine) are evaluated too: a plain ORF search
    must stop at that UGA, so they are expected misses and are reported
    separately.

    Matching is by the stop codon: an annotated CDS is *found* when an
    ORF on the same strand ends at the same base (the standard criterion
    for gene-finder sensitivity, because the stop is unambiguous while
    the start is a choice among in-frame starts). For found CDS the start
    is compared as well: identical, ORF start upstream (longer ORF) or
    downstream.

    Precision is not reported: an ORF search is not a gene predictor, and
    most long ORFs on the antisense strand of real genes are not genes.
    The share of ORFs >= 300 nt that coincide with an annotated stop is
    given as context only.

Usage::

    python scripts/benchmarks/orf_refseq.py [--json out.json]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from genome_scale import CACHE, GENOME_ACCESSION, genome_path  # noqa: E402

GFF_SHA256 = "f0c966bed5c950e6325fe1f885a1fbacee9144c30f5cd43feedf446768935617"
GFF_URL = ("https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
           f"?db=nuccore&report=gff3&id={GENOME_ACCESSION}&retmode=text")
MIN_LENGTH = 90


def gff_path() -> Path:
    path = CACHE / f"{GENOME_ACCESSION}.gff3"
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(GFF_URL, headers={"User-Agent": "synapse-benchmark/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != GFF_SHA256:
        raise SystemExit(f"ERROR: {path} has SHA-256 {digest}, expected {GFF_SHA256}.")
    return path


def annotated_cds() -> Tuple[List[dict], Counter]:
    """Return the evaluable CDS (1-based inclusive coordinates) and the excluded counts."""
    rows: Dict[str, List[List[str]]] = {}
    for line in gff_path().read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) == 9 and fields[2] == "CDS":
            identifier = fields[8].split(";", 1)[0]
            rows.setdefault(identifier, []).append(fields)
    excluded: Counter = Counter()
    cds = []
    for identifier, segments in rows.items():
        attributes = segments[0][8]
        if "pseudo=true" in attributes:
            excluded["pseudogene"] += 1
            continue
        if "exception=ribosomal slippage" in attributes:
            excluded["ribosomal slippage"] += 1
            continue
        if len(segments) > 1:
            excluded["multi-segment"] += 1
            continue
        start, end, strand = int(segments[0][3]), int(segments[0][4]), segments[0][6]
        gene = next((part[5:] for part in attributes.split(";") if part.startswith("gene=")), identifier)
        cds.append({"id": identifier, "gene": gene, "start": start, "end": end, "strand": strand,
                    "length": end - start + 1, "selenocysteine": "aa:Sec" in attributes})
    return cds, excluded


def synapse_orfs(genome: str) -> List[dict]:
    from core.api import find_orfs

    result = find_orfs(genome, min_length=MIN_LENGTH, start_codon_policy="atg_and_alternative",
                       nested=False, genetic_code=11, include_partial=False, top_n=10**9)
    return result.data["orfs"]


def evaluate(genome: str) -> dict:
    cds, excluded = annotated_cds()
    started = time.perf_counter()
    orfs = synapse_orfs(genome)
    runtime = time.perf_counter() - started
    # Synapse records: 0-based start of the first base, 0-based last base.
    by_stop = {}
    for orf in orfs:
        if orf["strand"] == "+":
            by_stop[("+", orf["stop"] + 1)] = orf["start"] + 1
        else:
            by_stop[("-", orf["start"] + 1)] = orf["stop"] + 1
    outcome: Counter = Counter()
    missed = []
    for entry in cds:
        # The stop of a plus-strand CDS is its end; of a minus-strand CDS its start.
        stop_key = (entry["strand"], entry["end"] if entry["strand"] == "+" else entry["start"])
        annotated_start = entry["start"] if entry["strand"] == "+" else entry["end"]
        group = "selenoprotein" if entry["selenocysteine"] else "regular"
        if stop_key not in by_stop:
            outcome[f"{group}:missed"] += 1
            missed.append({k: entry[k] for k in ("gene", "start", "end", "strand", "length")})
            continue
        orf_start = by_stop[stop_key]
        upstream = (orf_start < annotated_start) if entry["strand"] == "+" else (orf_start > annotated_start)
        if orf_start == annotated_start:
            outcome[f"{group}:start identical"] += 1
        elif upstream:
            outcome[f"{group}:ORF start upstream"] += 1
        else:
            outcome[f"{group}:ORF start downstream"] += 1
    regular = [c for c in cds if not c["selenocysteine"]]
    found = sum(v for k, v in outcome.items() if k.startswith("regular:") and k != "regular:missed")
    annotated_stops = {(c["strand"], c["end"] if c["strand"] == "+" else c["start"]) for c in cds}
    long_orfs = [key for key, _ in by_stop.items()]
    long_with_stop = sum(1 for orf in orfs if orf["length"] >= 300 and (
        (("+", orf["stop"] + 1) if orf["strand"] == "+" else ("-", orf["start"] + 1)) in annotated_stops))
    long_total = sum(1 for orf in orfs if orf["length"] >= 300)
    return {
        "genome": GENOME_ACCESSION,
        "settings": {"genetic_code": 11, "start_codons": ["ATG", "GTG", "TTG"], "nested": False,
                     "min_length_nt": MIN_LENGTH, "partial": False},
        "annotated_cds_evaluated": len(cds),
        "excluded": dict(excluded),
        "orfs_reported": len(long_orfs),
        "orf_search_seconds": round(runtime, 2),
        "outcome": dict(sorted(outcome.items())),
        "sensitivity_regular": round(found / len(regular), 4),
        "missed_regular_shorter_than_min": sum(1 for m in missed if m["length"] < MIN_LENGTH),
        "orfs_ge_300nt": long_total,
        "orfs_ge_300nt_ending_at_annotated_stop": long_with_stop,
        "missed": missed,
    }


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--json")
    args = parser.parse_args(argv)
    from core.api import load_records

    genome = load_records(str(genome_path()))[0].sequence
    report = evaluate(genome)
    summary = {k: v for k, v in report.items() if k != "missed"}
    summary["missed_examples"] = report["missed"][:25]
    print(json.dumps(summary, indent=2))
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
