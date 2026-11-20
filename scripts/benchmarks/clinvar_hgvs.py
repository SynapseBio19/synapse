"""Build the ClinVar HGVS benchmark set (tests/fixtures/clinvar/).

ClinVar (NCBI) records each variant with its HGVS name on a RefSeq
transcript, e.g. ``NM_000518.5(HBB):c.20A>T (p.Glu7Val)``. This script
collects such names for a few genes with short coding sequences, keeps the
exonic coding changes (no intronic offsets, no UTR positions) and stores
them together with the coding sequence of the transcript version they
refer to. ``tests/test_hgvs_clinvar_benchmark.py`` then rebuilds every
observed coding sequence and checks that Synapse describes it with the
same HGVS c. name (and, for codon substitutions, the same protein change).

ClinVar and RefSeq records are produced by NCBI, a U.S. government
agency, and are in the public domain; the accession list and query date
are written into the fixture so the selection can be traced.

The set is a regression benchmark for notation, not a clinical
validation: Synapse does not interpret pathogenicity.

Usage (needs network; run rarely)::

    python scripts/benchmarks/clinvar_hgvs.py
"""

from __future__ import annotations

import datetime
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tests" / "fixtures" / "clinvar"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HEADERS = {"User-Agent": "synapse-benchmark/1.0"}

# Genes with coding sequences short enough for a full global alignment
# (Synapse recommends <= 2,000 nt) and many ClinVar submissions.
GENES = ["HBB", "KRAS", "VHL", "TP53", "PTEN", "MECP2"]
PER_GENE = 250
SEARCH_LIMIT = 4000

# Exonic coding change without intronic offset (+/-) or UTR (*, -).
SIMPLE_C = re.compile(r"^c\.(\d+)(?:_(\d+))?(?:([ACGT])>([ACGT])|del|dup|ins([ACGT]+)|delins([ACGT]+))$")
TITLE = re.compile(r"^(N[MR]_\d+\.\d+)\(([^)]+)\):(c\.[^ ]+)(?: \((p\.[^)]+)\))?$")


def _get(endpoint: str, **params: str) -> bytes:
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=60) as response:
                payload = response.read()
            time.sleep(0.4)  # stay below NCBI's three requests per second
            return payload
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 + attempt * 3)
    raise RuntimeError("unreachable")


def clinvar_titles(gene: str) -> List[dict]:
    search = json.loads(_get("esearch.fcgi", db="clinvar", term=f"{gene}[gene]",
                             retmax=str(SEARCH_LIMIT), retmode="json"))
    ids = sorted(search["esearchresult"]["idlist"], key=int)
    titles: Dict[str, str] = {}
    for offset in range(0, len(ids), 200):
        titles.update(_summaries(ids[offset:offset + 200], gene))
    records = []
    for uid in sorted(titles, key=int):
        match = TITLE.match(titles[uid])
        if match and match.group(2) == gene:
            records.append({"uid": int(uid), "transcript": match.group(1),
                            "c": match.group(3), "p": match.group(4)})
    return records


def _summaries(batch: List[str], gene: str) -> Dict[str, str]:
    """UID -> title; a batch too large for NCBI's JSON conversion is split in half."""
    summary: dict = {}
    for attempt in range(5):
        summary = json.loads(_get("esummary.fcgi", db="clinvar", id=",".join(batch), retmode="json"))
        if "result" in summary:
            return {uid: summary["result"][uid].get("title", "")
                    for uid in summary["result"].get("uids", [])}
        message = str(summary.get("eutilsresult", {}).get("ERROR", ""))
        if "cannot be transformed to JSON" in message:
            if len(batch) == 1:
                return {}          # one oversized record: skipped, it is not a short coding change
            middle = len(batch) // 2
            return {**_summaries(batch[:middle], gene), **_summaries(batch[middle:], gene)}
        # NCBI answers an overload with {"error": ...}; wait and ask again.
        time.sleep(5 + 5 * attempt)
    raise SystemExit(f"ERROR: ClinVar esummary kept failing for {gene}: {summary}")


def coding_sequence(transcript: str) -> str:
    text = _get("efetch.fcgi", db="nuccore", id=transcript, rettype="fasta_cds_na", retmode="text").decode()
    entries = text.strip().split("\n>")
    if len(entries) != 1:
        raise SystemExit(f"ERROR: {transcript} has {len(entries)} CDS entries.")
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    selected: List[dict] = []
    transcripts: Dict[str, dict] = {}
    for gene in GENES:
        records = clinvar_titles(gene)
        transcript, _ = Counter(r["transcript"] for r in records).most_common(1)[0]
        cds = coding_sequence(transcript)
        if len(cds) % 3 or cds[:3] != "ATG":
            raise SystemExit(f"ERROR: unexpected CDS for {transcript}.")
        transcripts[transcript] = {"gene": gene, "cds_length": len(cds), "cds": cds}
        kept = 0
        for record in records:
            match = SIMPLE_C.match(record["c"])
            if record["transcript"] != transcript or not match:
                continue
            last = int(match.group(2) or match.group(1))
            if last > len(cds):
                continue
            selected.append({"uid": record["uid"], "gene": gene, "transcript": transcript,
                             "c": record["c"], "p": record["p"]})
            kept += 1
            if kept == PER_GENE:
                break
        print(f"{gene}: {transcript}, {len(records)} titles, {kept} kept", flush=True)
    fixture = {
        "source": "NCBI ClinVar (E-utilities esearch/esummary) and RefSeq (efetch fasta_cds_na)",
        "licence": "public domain (U.S. government work, NCBI)",
        "retrieved": datetime.date.today().isoformat(),
        "query": "<gene>[gene], lowest ClinVar UIDs first",
        "selection": "exonic coding changes on the gene's most frequent transcript; "
                     f"at most {PER_GENE} per gene",
        "transcripts": transcripts,
        "variants": selected,
    }
    path = OUT / "hgvs_benchmark.json"
    path.write_text(json.dumps(fixture, indent=1) + "\n", encoding="utf-8")
    print(f"{len(selected)} variants written to {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
