"""Genome-scale benchmark: every analysis on a complete bacterial genome.

The reference genome is *Escherichia coli* K-12 substr. MG1655
(RefSeq NC_000913.3, 4,641,652 bp). It is downloaded once from NCBI
E-utilities into ``.cache/benchmarks/`` (ignored by Git) and verified by
its SHA-256 checksum; nothing is uploaded.

Each case runs in its own child process with a time limit, so a case that
is slow or runs out of memory cannot hide the results of the others. The
report gives wall-clock time, CPU time and the peak resident memory of
the child process.

Usage (from the project folder)::

    python scripts/benchmarks/genome_scale.py              # all cases
    python scripts/benchmarks/genome_scale.py --case orf   # one case
    python scripts/benchmarks/genome_scale.py --json out.json

The numbers depend on the machine; record them together with the
interpreter version and the CPU (docs/BENCHMARKS.md).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Callable, Dict, List

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / ".cache" / "benchmarks"

GENOME_ACCESSION = "NC_000913.3"
GENOME_LENGTH = 4_641_652
GENOME_SHA256 = "6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7"
GENOME_URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
              f"?db=nuccore&id={GENOME_ACCESSION}&rettype=fasta&retmode=text")

DEFAULT_TIMEOUT = 900


def genome_path() -> Path:
    """Return the cached genome FASTA, downloading and verifying it if needed."""
    path = CACHE / f"{GENOME_ACCESSION}.fasta"
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(GENOME_URL, headers={"User-Agent": "synapse-benchmark/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
        path.write_bytes(payload)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != GENOME_SHA256:
        raise SystemExit(f"ERROR: {path} has SHA-256 {digest}, expected {GENOME_SHA256}. "
                         "Delete the file and run again.")
    return path


def _genome() -> str:
    sys.path.insert(0, str(ROOT))
    from core.api import load_records

    records = load_records(str(genome_path()))
    sequence = records[0].sequence
    assert len(sequence) == GENOME_LENGTH, len(sequence)
    return sequence


# ----------------------------------------------------------------------
# Cases. Each returns a short, comparable summary of its result.
# ----------------------------------------------------------------------
def case_load() -> str:
    return f"{len(_genome())} bp"


def case_overview() -> str:
    from core.api import sequence_overview

    result = sequence_overview(_genome())
    return f"GC {result.data['gc_content']:.2f}%" if isinstance(result.data, dict) else "ok"


def case_orf() -> str:
    from core.api import find_orfs

    result = find_orfs(_genome())
    data = result.data
    return f"{data.get('total_orfs', len(data.get('orfs', [])))} ORFs"


def case_rna() -> str:
    from core.api import analyze_rna

    analyze_rna(_genome())
    return "ok"


def case_restriction() -> str:
    from core.api import restriction_digest
    from data.restriction_enzymes import RESTRICTION_ENZYMES

    # Every enzyme of the database, including the frequent 4-bp cutters
    # (AluI, HaeIII, MspI, TaqI ...): the worst case for the fragment lists.
    result = restriction_digest(_genome(), sorted(RESTRICTION_ENZYMES))
    data = result.data
    sites = data.get("total_sites") if isinstance(data, dict) else None
    return f"{sites} sites" if sites is not None else "ok"


def case_primer() -> str:
    from core.api import design_primers

    result = design_primers(_genome(), 2_000_000, 2_000_600)
    return f"{len(result.data.ranked_pairs)} pairs"


def case_protein() -> str:
    from core.api import analyze_protein_from_cds

    # thrL..thrA region: thrA (b0002) CDS 337..2799 on the plus strand.
    analyze_protein_from_cds(_genome()[336:2799])
    return "ok"


def case_compare() -> str:
    from core.api import compare_sequences

    genome = _genome()
    reference = genome[:2_000_000]
    variant = reference[:1_000_000] + ("A" if reference[1_000_000] != "A" else "C") + reference[1_000_001:]
    result = compare_sequences(reference, variant)
    return f"{len(getattr(result.data, 'mutations', []))} differences"


CASES: Dict[str, Callable[[], str]] = {
    "load": case_load,
    "overview": case_overview,
    "orf": case_orf,
    "rna": case_rna,
    "restriction": case_restriction,
    "primer": case_primer,
    "protein": case_protein,
    "compare": case_compare,
}


# ----------------------------------------------------------------------
# Child-process runner
# ----------------------------------------------------------------------
def _peak_memory_bytes() -> int:
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        kernel32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        kernel32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
        counters = Counters()
        counters.cb = ctypes.sizeof(Counters)
        if not kernel32.K32GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            raise OSError(ctypes.get_last_error(), "K32GetProcessMemoryInfo failed")
        return int(counters.PeakWorkingSetSize)
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _child(name: str) -> None:
    sys.path.insert(0, str(ROOT))
    _genome()  # load outside the timed region, so "load" is measured on its own
    wall, cpu = time.perf_counter(), time.process_time()
    summary = CASES[name]()
    record = {"case": name, "wall_s": round(time.perf_counter() - wall, 2),
              "cpu_s": round(time.process_time() - cpu, 2),
              "peak_mib": round(_peak_memory_bytes() / 2**20, 1), "summary": summary}
    print("RESULT " + json.dumps(record))


def run_case(name: str, timeout: int) -> dict:
    command = [sys.executable, __file__, "--child", name]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=timeout, cwd=ROOT)
    except subprocess.TimeoutExpired:
        return {"case": name, "status": "timeout", "timeout_s": timeout}
    for line in completed.stdout.splitlines():
        if line.startswith("RESULT "):
            return {"status": "ok", **json.loads(line[7:])}
    tail = (completed.stderr or completed.stdout).strip().splitlines()[-3:]
    return {"case": name, "status": "error", "detail": " | ".join(tail)}


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--json")
    parser.add_argument("--child", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.child:
        _child(args.child)
        return 0
    genome_path()
    results = []
    for name in args.case or list(CASES):
        result = run_case(name, args.timeout)
        results.append(result)
        print(json.dumps(result), flush=True)
    meta = {"python": platform.python_version(), "machine": platform.machine(),
            "processor": platform.processor(), "system": platform.platform(),
            "genome": GENOME_ACCESSION, "os_cpu_count": os.cpu_count()}
    if args.json:
        Path(args.json).write_text(json.dumps({"meta": meta, "results": results}, indent=2), encoding="utf-8")
    return 0 if all(r["status"] == "ok" for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
