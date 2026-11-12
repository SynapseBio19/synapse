"""Stage F, recommendation 5: every menu of the interface, rendered and driven.

Each scenario runs a real dialog with scripted keystrokes and records every
``style.menu`` call it draws (call site, entries, ``[0]`` label, the text on
screen). The checks cover entries and numbering, navigation with ``0``,
invalid input, cancelling, and the error states and result views behind the
entries. A last check compares the recorded call sites with every
``style.menu`` call in the sixteen interface files, so a new menu without a
rendered test fails the suite.

No network: the NCBI client is replaced by a parser of a fixed BLAST XML.
"""

from __future__ import annotations

import ast
import builtins
import contextlib
import functools
import io
import logging
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests._support import PROJECT_ROOT, deterministic_dna
from tests.test_terminal_ui import INTERFACE_SOURCES
from ui import style
from ui import blast_menu, cloning_menu, crispr_menu, mutation_menu, orf_menu, primer_menu, restriction_menu, variant_analysis_menu

DNA = deterministic_dna(900)
# Vector and insert with EcoRI (GAATTC) and BamHI (GGATCC) sites for the cloning dialogs.
VECTOR = deterministic_dna(400) + "GAATTC" + deterministic_dna(300)[::-1] + "GGATCC" + deterministic_dna(200)
INSERT = "GAATTC" + deterministic_dna(240) + "GGATCC"
ZERO_LABELS = {style.BACK_LABEL, style.CANCEL_LABEL, style.EXIT_LABEL, "Cancel (choose another sequence)"}


class _Exhausted(Exception):
    """The script ran out of answers before the dialog returned."""


class Session:
    """Run one dialog with scripted answers; record its output and every menu it draws."""

    def __init__(self, answers):
        self.pending = list(answers)
        self.buffer = io.StringIO()
        self.menus = []
        self.screen = ""
        self.result = None
        self.exit = None

    def _input(self, prompt=""):
        self.buffer.write(prompt)
        if not self.pending:
            raise _Exhausted(prompt)
        answer = self.pending.pop(0)
        self.buffer.write(answer + "\n")
        return answer

    def _menu(self, groups, *, zero=style.BACK_LABEL, _original=style.menu):
        caller = sys._getframe(1)
        start = self.buffer.tell()
        _original(groups, zero=zero)
        self.menus.append({
            "site": (Path(caller.f_code.co_filename).resolve().relative_to(PROJECT_ROOT).as_posix(),
                     caller.f_lineno),
            "entries": [(str(key), label) for _title, entries in groups for key, label in entries],
            "titles": [title for title, _entries in groups],
            "zero": zero,
            "screen": self.screen,
            "text": self.buffer.getvalue()[start:],
        })

    def _header(self, title, *args, _original=style.header, **kwargs):
        self.screen = str(title).upper()
        _original(title, *args, **kwargs)

    def run(self, call):
        logging.disable(logging.CRITICAL)
        try:
            with mock.patch.object(builtins, "input", self._input), \
                    mock.patch.object(style, "menu", self._menu), \
                    mock.patch.object(style, "header", self._header), \
                    contextlib.redirect_stdout(self.buffer):
                try:
                    self.result = call()
                except SystemExit as exc:
                    self.exit = exc
        finally:
            logging.disable(logging.NOTSET)
        self.output = self.buffer.getvalue()
        return self

    def menus_on(self, screen):
        """The menus drawn under the screen title `screen`."""
        return [m for m in self.menus if m["screen"] == screen]


def menu_call_sites():
    """Every style.menu(...) call in the interface files, as (path, first line)."""
    sites = set()
    for path in INTERFACE_SOURCES:
        rel = path.relative_to(PROJECT_ROOT).as_posix()
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                    and node.func.attr == "menu" and isinstance(node.func.value, ast.Name) \
                    and node.func.value.id == "style":
                sites.add((rel, node.lineno))
    return sites


# ----------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------

def _fasta(directory):
    path = os.path.join(directory, "two.fasta")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(">rec1\n" + DNA[:300] + "\n>rec2\n" + DNA[300:600] + "\n")
    return path


def _fake_ncbi_run(self, query, parameters):
    from analyses import blast
    from tests.test_p0_regressions import _blast_xml
    return blast.parse_blast_xml(_blast_xml(1, 32, 32, 1e-5, "A" * 32, "A" * 32), query=query, rid="T",
                                 parameters=parameters, submitted_at="2026-01-01T00:00:00+00:00")


def _local_databases():
    from analyses import blast
    return [blast.LocalBlastDatabase(name="refseq_rna", database_type="nucleotide", path="x",
                                     installed=True, structurally_valid=True,
                                     verification_status="VERIFIED", size_bytes=3 * 2 ** 30,
                                     sequence_count=12345),
            blast.LocalBlastDatabase(name="custom", database_type="nucleotide", path="y",
                                     installed=True, structurally_valid=True,
                                     verification_status="USER_SUPPLIED", size_bytes=2 ** 20)]


def _orf_settings():
    return {"min_length": 75, "start_codon_policy": "atg_only", "nested": False,
            "ranking": "length", "include_partial": True, "genetic_code": 1}


@functools.lru_cache(maxsize=None)
def scenarios():
    """Run every scenario once; the tests below read the recorded sessions."""
    from analyses import blast, cloning, crispr, mutation, orf, primer, restriction
    from core.models import SequenceRecord
    import main
    from ui import menu as ui_menu

    work = tempfile.mkdtemp(prefix="synapse_menus_")
    fasta = _fasta(work)
    orfs = orf.analyze_orfs(DNA, min_length=60)["orfs"][:2]
    records = [{"header": "ref", "sequence": DNA[:120]}, {"header": "var", "sequence": DNA[:60] + "T" + DNA[61:120]}]
    manager_root = os.path.join(work, "blastdb")
    runs = {}

    def run(name, answers, call, patches=()):
        with contextlib.ExitStack() as stack:
            for target, attribute, value in patches:
                stack.enter_context(mock.patch.object(target, attribute, value))
            try:
                runs[name] = Session(answers).run(call)
            except _Exhausted as exc:  # pragma: no cover - a broken script is a test failure
                raise AssertionError(f"scenario {name}: script exhausted at prompt {exc}") from exc

    no_network = [(blast.NCBIBlastClient, "run", _fake_ncbi_run)]
    with_databases = [(blast.LocalBlastDatabaseManager, "list_databases", lambda self: _local_databases())]
    without_databases = [(blast.LocalBlastDatabaseManager, "list_databases", lambda self: [])]

    # --- BLAST database manager: menu, type menu, NCBI catalog, status and tool views
    run("blast_manager", ["9", "4", "x.fa", "db", "9", "4", "x.fa", "db", "0", "2", "x", "2", "0",
                          "1", "7", "0"],
        lambda: blast_menu.blast_database_manager_menu(blast.LocalBlastDatabaseManager(manager_root)))
    # --- BLAST search: query source, FASTA records, engine, program, remote database, hits
    run("blast_search", ["1", "9", "1", "0",
                         "1", "2", fasta, "9",
                         "1", "2", fasta, "0",
                         "1", "2", fasta, "1", "9",
                         "1", "2", fasta, "1", "2", "9",
                         "1", "2", fasta, "1", "2", "1", "0",
                         "1", "2", fasta, "1", "2", "1", "1", "", "", "y",
                         "6", "0", "6", "1", "7", "0"],
        lambda: blast_menu.blast_menu(DNA, "demo"), no_network)
    # --- BLAST local engine: the local database list, invalid choice, cancel, no database installed
    run("blast_local", ["1", "2", fasta, "1", "1", "1", "9",
                        "1", "2", fasta, "1", "1", "1", "0", "0"],
        lambda: blast_menu.blast_menu(DNA, "demo"), with_databases)
    run("blast_local_none", ["1", "2", fasta, "1", "1", "1", "0"],
        lambda: blast_menu.blast_menu(DNA, "demo"), without_databases)
    # --- BLAST query type for a sequence whose letters fit DNA and protein
    for answer in ("1", "2", "0", "9"):
        run(f"blast_query_type_{answer}", [answer], lambda: blast_menu._query_type_for_menu("ATGCRYATGC", None))

    # --- Restriction: analysis menu, enzyme input errors, topology, every view
    run("restriction", ["9",
                        "1", "FooI", "EcoRI", "", "",
                        "2", "EcoRI", "", "",
                        "3", "EcoRI,BamHI", "0",
                        "3", "EcoRI,BamHI", "2", "",
                        "3", "BsaI,EcoRI", "", "",
                        "4", "", "EcoRI", "", "",
                        "0"],
        lambda: restriction_menu.run_restriction_analysis_ui(VECTOR, header="vec"))
    run("restriction_topology", ["x", ""], restriction_menu._synapse_ui_select_topology)
    run("restriction_topology_circular", ["2"], restriction_menu._synapse_ui_select_topology)
    run("restriction_topology_back", ["0"], restriction_menu._synapse_ui_select_topology)

    # --- Cloning: enzyme list with and without [0], single-enzyme and two-enzyme dialogs
    run("cloning_enzyme", ["x", "99", "0"], lambda: cloning_menu._cloning_ui_select_enzyme("Selection"))
    run("cloning_enzyme_required", ["0", "2"],
        lambda: cloning_menu._cloning_ui_select_enzyme("Selection", allow_back=False))
    run("cloning_single", ["9", "1", "1", ""], lambda: cloning_menu._cloning_ui_single_enzyme(VECTOR, INSERT))
    run("cloning_single_back", ["0"], lambda: cloning_menu._cloning_ui_single_enzyme(VECTOR, INSERT))
    run("cloning_two", ["9", "1", "1", "1", "", "1", "2", "1", "2", "9", "1", ""],
        lambda: cloning_menu._cloning_ui_two_enzyme_strategy(VECTOR, INSERT))
    run("cloning_two_cancel", ["1", "1", "2", "1", "2", "0"],
        lambda: cloning_menu._cloning_ui_two_enzyme_strategy(VECTOR, INSERT))
    run("cloning_planner", ["9", "1", ""], lambda: cloning_menu._cloning_ui_automatic_planner(VECTOR, INSERT))
    run("cloning_planner_back", ["0"], lambda: cloning_menu._cloning_ui_automatic_planner(VECTOR, INSERT))
    candidates = [{"header": "vector", "sequence": VECTOR}, {"header": "insert", "sequence": INSERT},
                  {"header": "empty", "sequence": ""}]
    run("cloning_insert", ["9", "1"], lambda: cloning_menu._cloning_ui_select_insert(candidates, vector_index=0))
    run("cloning_insert_back", ["0"], lambda: cloning_menu._cloning_ui_select_insert(candidates, vector_index=0))
    # the cloning menu itself, reaching the single-enzyme analysis through it
    run("cloning_menu", ["9", "2", "1", "1", "", "0"], lambda: cloning_menu.run_cloning_analysis_ui(VECTOR, INSERT))

    # --- CRISPR: every sub-menu, the two option lists, the batch records
    run("crispr", ["5", "11",                                   # notes before any result
                   "1", "9", "1", "manual", DNA[:200], "", "2", "1", "120", "0",
                   "2", "2", "", "1", fasta, "2", "9", "2", "1", "0",
                   "3", "9", "1", "9", "2", "2", "0",
                   "4", "", "",
                   "5", "9", "2", "1", "1", "0",
                   "11", "9", "1", "2", "guides.txt", "2", os.path.join(work, "guides.json"),
                   "3", os.path.join(work, "batch.csv"), "0", "0"],
        lambda: crispr_menu.crispr_menu(DNA, "demo"))
    run("crispr_editor", ["x", "0"], crispr_menu._choose_editor)
    run("crispr_system_cancel", ["0"], crispr_menu._choose_system)

    # --- Mutation: sub-menu, the sequence prompt without [0], manual entry
    run("mutation", ["9", "1", "9", "1", "2", "", "0"],
        lambda: mutation_menu.run_mutation_analysis_ui(sequences=records))
    # [0] at the reference prompt, and at the variant prompt after a reference was chosen
    run("mutation_cancel", ["1", "0", "1", "1", "0", "0"],
        lambda: mutation_menu.run_mutation_analysis_ui(sequences=records))

    # --- ORF: selector and settings
    run("orf_select", ["9", "1"], lambda: orf_menu.select_orf(orfs))
    run("orf_select_back", ["0"], lambda: orf_menu.select_orf(orfs))
    run("orf_settings", ["9", "2", "3", "4", "5", "1", "x", "1", "100", "6", "11", "0"],
        lambda: orf_menu._orf_settings_menu(_orf_settings()))
    run("orf_settings_restore", ["8", "7", "0"], lambda: orf_menu._orf_settings_menu(_orf_settings()))

    # --- Primer: Tm model
    run("primer_tm", ["9", "2"], primer_menu._primer_ui_select_tm_model)
    run("primer_tm_cancel", ["0"], primer_menu._primer_ui_select_tm_model)
    # --- Primer: design around a restriction site (workflow 10.2)
    run("primer_site", ["ecori", "9", "1", "", "1", "", "", "", ""],
        lambda: primer_menu._primer_ui_restriction_site_design(VECTOR, None))
    run("primer_site_cancel", ["EcoRI", "0"], lambda: primer_menu._primer_ui_restriction_site_design(VECTOR, None))
    run("primer_site_unknown", ["NoSuchI"], lambda: primer_menu._primer_ui_restriction_site_design(VECTOR, None))
    run("primer_site_none", ["SmaI"], lambda: primer_menu._primer_ui_restriction_site_design(VECTOR, None))

    # --- Session: type declaration, protein record, sequence selector
    ambiguous = SequenceRecord(identifier="amb", sequence="ATGCRYATGC")
    run("main_type", ["9", "1"], lambda: main.resolve_sequence_type(ambiguous))
    run("main_type_cancel", ["0"], lambda: main.resolve_sequence_type(ambiguous))
    protein = SequenceRecord(identifier="pep", sequence="MKTAYIAKQRQISFVKSHFSRQ", alphabet="PROTEIN")
    run("main_protein", ["9", "", "1", "0", "S"], lambda: main.run_protein_record(protein, 1))
    run("main_protein_exit", ["0"], lambda: main.run_protein_record(protein, 1))
    loaded = [{"header": "a", "sequence": DNA[:90]}, {"header": "b", "sequence": DNA[:60], "topology": "circular"}]
    run("selector", ["9", "x", "2"], lambda: ui_menu.select_gene(loaded))
    run("selector_file", ["F"], lambda: ui_menu.select_gene(loaded))
    return runs


def _menu_block_problems(menu):
    """Layout and numbering problems of one recorded menu."""
    problems = []
    keys = [key for key, _label in menu["entries"]]
    numbers = [int(key) for key in keys if key.isdigit()]
    if len(set(keys)) != len(keys):
        problems.append(f"duplicate keys {keys}")
    if numbers and numbers != list(range(1, len(numbers) + 1)):
        problems.append(f"numbers not 1..n: {numbers}")
    if menu["zero"] is not None and menu["zero"] not in ZERO_LABELS:
        problems.append(f"[0] means {menu['zero']!r}")
    lines = menu["text"].splitlines()
    if lines[-1] != style.minor_rule():
        problems.append("does not end with the minor rule")
    shown = [line for line in lines if line.startswith("[")]
    expected = [f"[{key}] {label}" for key, label in menu["entries"]]
    if menu["zero"] is not None:
        expected.append(f"[0] {menu['zero']}")
    if shown != expected:
        problems.append(f"entries on screen {shown} != {expected}")
    if menu["zero"] is not None and lines[-3] != "":
        problems.append("no blank line before [0]")
    return problems


class EveryMenuIsRenderedTests(unittest.TestCase):
    def test_every_menu_call_site_is_rendered_by_a_scenario(self):
        drawn = {menu["site"] for session in scenarios().values() for menu in session.menus}
        drawn |= _sites_rendered_by_the_screen_tests()
        missing = sorted(menu_call_sites() - drawn)
        self.assertEqual(missing, [], "menus without a rendered test")

    def test_every_rendered_menu_has_the_shared_layout(self):
        for name, session in scenarios().items():
            for menu in session.menus:
                with self.subTest(scenario=name, site=menu["site"]):
                    self.assertEqual(_menu_block_problems(menu), [])


@functools.lru_cache(maxsize=None)
def _sites_rendered_by_the_screen_tests():
    """The module menus that RenderedScreenTests already draw (first screens)."""
    from tests.test_terminal_ui import LEAD_IN, SHARED_UI_MODULES, render_first_screen
    sites = set()
    original = style.menu

    def record(groups, *, zero=style.BACK_LABEL):
        caller = sys._getframe(1)
        sites.add((Path(caller.f_code.co_filename).resolve().relative_to(PROJECT_ROOT).as_posix(),
                   caller.f_lineno))
        original(groups, zero=zero)

    with mock.patch.object(style, "menu", record):
        for module in SHARED_UI_MODULES:
            render_first_screen(module, LEAD_IN.get(module, ()) + ("0",))
    return sites


class BlastMenuTests(unittest.TestCase):
    def test_database_manager(self):
        s = scenarios()["blast_manager"]
        main_menu = s.menus[0]
        self.assertEqual([k for k, _ in main_menu["entries"]], [str(n) for n in range(1, 9)])
        self.assertEqual(main_menu["titles"], ["Database", "Tools"])
        type_menus = [m for m in s.menus if m["titles"] == ["Database Type"]]
        self.assertEqual(len(type_menus), 2)
        self.assertEqual(type_menus[0]["entries"], [("1", "Nucleotide"), ("2", "Protein")])
        self.assertEqual(type_menus[0]["zero"], style.CANCEL_LABEL)
        catalog = [m for m in s.menus if m["titles"] == ["NCBI databases available through update_blastdb.pl"]]
        self.assertEqual(len(catalog), 2)
        self.assertEqual(len(catalog[0]["entries"]), 6)
        # unknown key, invalid type, invalid catalog entry: three replies, each back to the manager
        self.assertEqual(s.output.count(style.INVALID_SELECTION), 3)
        self.assertIn("LOCAL BLAST DATABASE STATUS", s.output)     # [1] status view
        self.assertIn("BLAST+ TOOL STATUS", s.output)              # [7] tool view
        self.assertIsNone(s.result)

    def test_search_dialog_lists(self):
        s = scenarios()["blast_search"]
        by_title = {}
        for m in s.menus:
            by_title.setdefault(tuple(m["titles"]), []).append(m)
        source = by_title[("Query Source",)][0]
        self.assertEqual(source["entries"], [("1", "Selected Synapse Sequence"), ("2", "Load FASTA File")])
        records = by_title[("FASTA records",)][0]
        self.assertEqual(records["entries"], [("1", "rec1 | 300 dna"), ("2", "rec2 | 300 dna")])
        self.assertEqual(records["zero"], style.CANCEL_LABEL)
        self.assertEqual(by_title[("Search Engine",)][0]["entries"],
                         [("1", "Local BLAST+"), ("2", "NCBI BLAST")])
        programs = by_title[("BLAST program",)][0]
        self.assertEqual([label.split(" - ")[0] for _, label in programs["entries"]],
                         ["BLASTN", "BLASTX", "TBLASTX"])
        remote = by_title[("NCBI remote database",)][0]
        self.assertEqual([label.split(" - ")[0] for _, label in remote["entries"]],
                         ["nt", "refseq_rna", "refseq_genomic"])
        hits = by_title[("Ranked hits",)]
        self.assertEqual(hits[0]["entries"], [("1", "X1 | 100.00% identity | 10.67% coverage")])
        self.assertEqual(hits[0]["zero"], style.CANCEL_LABEL)

    def test_search_dialog_invalid_input_and_cancel(self):
        out = scenarios()["blast_search"].output
        # invalid source, FASTA record, engine and program, each answered and back to the menu
        self.assertEqual(out.count(style.INVALID_SELECTION), 4)
        self.assertIn("Database: 0", out)                            # cancel at the remote database
        self.assertIn("OK: BLAST search completed.", out)            # the mocked remote run
        self.assertIn("BLAST RESULT SUMMARY", out)
        self.assertIn("Hit number: 1", out)
        self.assertIn("Subject         1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  32", out)   # alignment view
        self.assertIn("1. X1: ", out)                                # evidence profile, a result list
        self.assertEqual(scenarios()["blast_search"].result.hits[0].accession, "X1")

    def test_local_database_list(self):
        s = scenarios()["blast_local"]
        local = [m for m in s.menus if m["titles"] == ["Local database"]]
        self.assertEqual(len(local), 2)
        self.assertEqual(local[0]["entries"], [("1", "refseq_rna | 3.00 GB | VERIFIED | 12,345 seq"),
                                               ("2", "custom | 0.00 GB | USER_SUPPLIED")])
        self.assertEqual(local[0]["zero"], style.CANCEL_LABEL)
        self.assertIn(style.INVALID_SELECTION, s.output)
        none = scenarios()["blast_local_none"].output
        self.assertIn("NOTE: No compatible local BLAST database is installed.", none)
        self.assertIn("NOTE: Use [4] BLAST Database Manager to download or create one.", none)


class BlastQueryTypeTests(unittest.TestCase):
    def test_undetermined_query_type_is_asked(self):
        s = scenarios()["blast_query_type_1"]
        self.assertEqual(s.menus[0]["entries"], [("1", "DNA"), ("2", "Protein")])
        self.assertEqual(s.menus[0]["zero"], style.CANCEL_LABEL)
        self.assertIn("cannot be determined", s.output)
        self.assertEqual(s.result, "DNA")
        self.assertEqual(scenarios()["blast_query_type_2"].result, "PROTEIN")
        # 0 cancels; an unknown key cancels too (the caller returns to the BLAST menu)
        self.assertIsNone(scenarios()["blast_query_type_0"].result)
        self.assertIsNone(scenarios()["blast_query_type_9"].result)


class RestrictionMenuTests(unittest.TestCase):
    def test_analysis_menu_and_views(self):
        s = scenarios()["restriction"]
        analysis = s.menus_on("RESTRICTION ANALYSIS")[0]
        self.assertEqual(analysis["entries"], [("1", "Restriction Overview"), ("2", "Restriction Map"),
                                               ("3", "Digest & Gel Analysis"), ("4", "Complete Restriction Report")])
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertIn("ERROR: Unknown enzyme(s): FooI", s.output)
        self.assertIn("ERROR: Please select at least 1 enzyme(s).", s.output)
        self.assertIn("23. BsaI       GGTCTC(1/5)", s.output)          # Type IIS offered in the digest
        self.assertIn("Enzymes:               BsaI, EcoRI", s.output)      # no BsaI site: EcoRI only
        for view in ("RESTRICTION OVERVIEW", "RESTRICTION MAP", "DIGEST & GEL ANALYSIS",
                     "COMPLETE RESTRICTION REPORT", "PREDICTED GEL BANDS", "FRAGMENTS"):
            with self.subTest(view=view):
                self.assertIn(f"\n{view}\n", s.output)
        self.assertIsNone(s.result)

    def test_topology_question(self):
        s = scenarios()["restriction_topology"]
        self.assertEqual(s.menus[0]["entries"], [("1", "Linear"), ("2", "Circular"),
                                                 ("Enter", "Keep current (Linear)")])
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertEqual(s.result, "linear")                           # Enter keeps the current one
        self.assertEqual(scenarios()["restriction_topology_circular"].result, "circular")
        self.assertIsNone(scenarios()["restriction_topology_back"].result)


class CloningMenuTests(unittest.TestCase):
    def test_cloning_menu(self):
        s = scenarios()["cloning_menu"]
        menu = s.menus_on("CLONING ANALYSIS")[0]
        self.assertEqual(menu["entries"], [("1", "Automatic Cloning Planner"), ("2", "Single-Enzyme Analysis"),
                                           ("3", "Two-Enzyme Strategy Analysis")])
        self.assertEqual(menu["titles"], ["Planning", "Analysis"])
        self.assertIn(style.INVALID_SELECTION, s.output)
        # [2] reaches the single-enzyme report, which used to raise KeyError: 'left'
        self.assertIn("Strategy:              Single-enzyme cloning", s.output)
        self.assertIn("Overall valid:         YES", s.output)

    def test_planner_and_insert_selection(self):
        planner = scenarios()["cloning_planner"]
        self.assertEqual(planner.menus[0]["entries"], [("1", "Find best strategy"), ("2", "Show top strategies")])
        self.assertIn(style.INVALID_SELECTION, planner.output)
        self.assertIn("RANKED STRATEGIES", planner.output)
        self.assertIsInstance(planner.result, dict)
        self.assertIsNone(scenarios()["cloning_planner_back"].result)
        insert = scenarios()["cloning_insert"]
        self.assertEqual(insert.menus[0]["entries"], [("1", "insert")])      # no vector, no empty record
        self.assertIn(style.INVALID_SELECTION, insert.output)
        self.assertEqual(insert.result[1], "insert")
        self.assertEqual(scenarios()["cloning_insert_back"].result, (None, None, None))

    def test_enzyme_list(self):
        s = scenarios()["cloning_enzyme"]
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        self.assertEqual([label for _, label in s.menus[0]["entries"]], list(RESTRICTION_ENZYMES))
        self.assertEqual(s.output.count(style.INVALID_SELECTION), 2)       # "x" and "99"
        self.assertIsNone(s.result)                                        # [0] Back
        required = scenarios()["cloning_enzyme_required"]
        self.assertIsNone(required.menus[0]["zero"])                       # no [0] offered
        self.assertIn(style.INVALID_SELECTION, required.output)           # so 0 is not a key
        self.assertEqual(required.result, list(RESTRICTION_ENZYMES)[1])

    def test_single_enzyme_dialog(self):
        s = scenarios()["cloning_single"]
        self.assertEqual(s.menus[0]["entries"], [("1", "Select enzyme")])
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertIsNotNone(s.result)                                     # the report ran
        self.assertIn("SINGLE-ENZYME ANALYSIS", s.output)
        self.assertIsNone(scenarios()["cloning_single_back"].result)

    def test_two_enzyme_dialog(self):
        s = scenarios()["cloning_two"]
        self.assertEqual(s.menus[0]["entries"], [("1", "Configure enzyme architecture")])
        run_menu = [m for m in s.menus if m["entries"] == [("1", "Run analysis")]]
        self.assertEqual(len(run_menu), 1)
        self.assertIn("ERROR: The two enzymes must be different.", s.output)
        self.assertIn("TWO-ENZYME STRATEGY", s.output)
        self.assertEqual(s.output.count(style.INVALID_SELECTION), 2)
        self.assertIsNotNone(s.result)
        self.assertIsNone(scenarios()["cloning_two_cancel"].result)


class CrisprMenuTests(unittest.TestCase):
    def test_sub_menus_and_their_entries(self):
        s = scenarios()["crispr"]
        expected = {
            "SEQUENCE & REGION": [("1", "Select DNA Sequence"), ("2", "Select ORF / Region")],
            "BATCH ANALYSIS": [("1", "Load FASTA / Batch Analysis"), ("2", "Select Record as Current Sequence")],
            "CRISPR SYSTEM": [("1", "Select CRISPR System"), ("2", "Show System Details")],
            "GUIDE ANALYSIS & RANKING": [("1", "Analyze Individual Guide"), ("2", "Rank Guides")],
            "REPORTS & EXPORT": [("1", "Full Guide Report"), ("2", "Export Guides"), ("3", "Export Batch Summary")],
        }
        for title, entries in expected.items():
            with self.subTest(menu=title):
                drawn = s.menus_on(title)
                self.assertTrue(drawn, f"{title} was not drawn")
                self.assertEqual(drawn[0]["entries"], entries)
                self.assertEqual(drawn[0]["zero"], style.BACK_LABEL)

    def test_notes_errors_and_results(self):
        out = scenarios()["crispr"].output
        self.assertIn("NOTE: Run [4] Target Discovery first.", out)          # [5] before discovery
        self.assertIn("NOTE: No analysis results available.", out)          # [11] before discovery
        self.assertIn("Selected sequence:     manual (200 bp)", out)
        self.assertIn("Selected Region 1-120 (120 bp).", out)
        self.assertIn("Current sequence:      rec1", out)
        self.assertIn("Selected system:       SpCas9-NG", out)
        self.assertIn("RANKED GUIDES (", out)                                # [5] → [2] Rank Guides
        self.assertIn("ERROR: Use .json, .csv or .fasta/.fa extension.", out)
        self.assertIn("OK: Exported: ", out)
        self.assertGreaterEqual(out.count(style.INVALID_SELECTION), 5)

    def test_option_lists(self):
        s = scenarios()["crispr"]
        systems = [m for m in s.menus if m["titles"] == ["Available CRISPR systems"]][0]
        self.assertEqual([label.split()[0] for _, label in systems["entries"]],
                         ["SpCas9", "SpCas9-NG", "SaCas9", "Cas12a"])
        self.assertEqual(systems["zero"], style.CANCEL_LABEL)
        records = [m for m in s.menus if m["titles"] == ["Available records"]]
        self.assertEqual(records[0]["entries"], [("1", "rec1"), ("2", "rec2")])
        editor = scenarios()["crispr_editor"]
        self.assertEqual([label.split()[0] for _, label in editor.menus[0]["entries"]],
                         ["CBE", "ABE", "CBE+", "ABE+"])
        self.assertIn(style.INVALID_SELECTION, editor.output)
        self.assertEqual(editor.result, "")                                  # [0] Cancel
        self.assertEqual(scenarios()["crispr_system_cancel"].result, "")


class MutationMenuTests(unittest.TestCase):
    def test_analysis_menu_and_sequence_prompt(self):
        s = scenarios()["mutation"]
        analysis = s.menus[0]
        self.assertEqual([k for k, _ in analysis["entries"]], ["1", "2", "3", "4", "5"])
        prompts = [m for m in s.menus if m["site"][1] != analysis["site"][1]]
        self.assertEqual(prompts[0]["entries"], [("1", "ref"), ("2", "var"), ("M", "Enter sequence manually")])
        self.assertEqual(prompts[0]["zero"], style.CANCEL_LABEL)
        self.assertGreaterEqual(s.output.count(style.INVALID_SELECTION), 2)
        self.assertIsNotNone(s.result)                           # the comparison ran
        self.assertIn("Variant selection: 2", s.output)

    def test_sequence_prompt_cancel_returns_to_the_menu(self):
        s = scenarios()["mutation_cancel"]
        analysis_menus = [m for m in s.menus if m["screen"] == "MUTATION ANALYSIS"
                          and m["entries"][0] == ("1", "Compare Sequences")]
        self.assertEqual(len(analysis_menus), 3)                 # drawn again after each cancel
        self.assertNotIn("ERROR", s.output)
        self.assertIn("Reference selection: 0", s.output)
        self.assertIn("Variant selection: 0", s.output)
        self.assertIsNone(s.result)                              # nothing was compared


class OrfMenuTests(unittest.TestCase):
    def test_selector(self):
        s = scenarios()["orf_select"]
        self.assertEqual([k for k, _ in s.menus[0]["entries"]], ["1", "2"])
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertIn("NOTE: Enter a candidate number or 0 to go back.", s.output)
        self.assertIs(s.result, scenarios()["orf_select"].result)
        self.assertEqual(s.result["rank"], 1)
        self.assertIsNone(scenarios()["orf_select_back"].result)

    def test_settings_entries_follow_the_values(self):
        s = scenarios()["orf_settings"]
        first, last = s.menus[0], s.menus[-1]
        self.assertEqual([k for k, _ in first["entries"]], [str(n) for n in range(1, 9)])
        self.assertIn("Start codons:        atg_only", dict(first["entries"])["2"])
        self.assertIn("atg_and_alternative", dict(last["entries"])["2"])
        self.assertIn("100 nt", dict(last["entries"])["1"])
        self.assertIn("ERROR: Please enter an integer >= 3.", s.output)
        self.assertEqual(s.result["genetic_code"], 11)
        # [5] cycles 3' end -> 5' and 3' ends -> excluded (one press here)
        self.assertEqual((s.result["nested"], s.result["ranking"], s.result["include_partial"],
                          s.result["include_5prime_partial"]), (True, "heuristic", True, True))
        restored = scenarios()["orf_settings_restore"].result
        self.assertEqual((restored["min_length"], restored["start_codon_policy"], restored["genetic_code"]),
                         (75, "atg_only", 1))


class PrimerMenuTests(unittest.TestCase):
    def test_tm_model(self):
        s = scenarios()["primer_tm"]
        self.assertEqual(s.menus[0]["entries"], [("1", "Nearest-neighbor"), ("2", "Wallace"), ("3", "Salt-adjusted")])
        self.assertEqual(s.menus[0]["zero"], style.CANCEL_LABEL)
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertEqual(s.result, "wallace")
        self.assertIsNone(scenarios()["primer_tm_cancel"].result)

    def test_design_around_a_restriction_site(self):
        s = scenarios()["primer_site"]
        self.assertEqual(s.menus[0]["entries"], [("1", "Position 401-406  GAATTC  (cut after 401)")])
        self.assertEqual(s.menus[0]["zero"], style.CANCEL_LABEL)
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertEqual((s.result.enzyme, s.result.site_start), ("EcoRI", 400))
        self.assertIn("RESTRICTION SITE AND DIAGNOSTIC DIGEST", s.output)
        self.assertIn("Fragments (bp)", s.output)
        self.assertIsNone(scenarios()["primer_site_cancel"].result)
        self.assertIn("ERROR: Unknown restriction enzyme: NoSuchI.", scenarios()["primer_site_unknown"].output)
        self.assertIn("ERROR: SmaI has no confirmed site in this sequence.", scenarios()["primer_site_none"].output)


class SessionMenuTests(unittest.TestCase):
    def test_sequence_type_declaration(self):
        s = scenarios()["main_type"]
        self.assertEqual(s.menus[0]["entries"], [("1", "DNA"), ("2", "RNA"), ("3", "Protein")])
        self.assertEqual(s.menus[0]["zero"], "Cancel (choose another sequence)")
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertEqual(s.result.alphabet, "DNA")
        self.assertIsNone(scenarios()["main_type_cancel"].result)

    def test_protein_record_menu(self):
        s = scenarios()["main_protein"]
        self.assertEqual(s.menus[0]["entries"], [("1", "Protein Analysis"), ("2", "BLAST Search"),
                                                 ("R", "Export Analysis Report (JSON / Markdown)"),
                                                 ("S", "Select Another Sequence")])
        self.assertEqual(s.menus[0]["zero"], style.EXIT_LABEL)
        self.assertIn(style.INVALID_SELECTION, s.output)
        self.assertIn("PROTEIN ANALYSIS", s.output)                       # [1] and back with 0
        self.assertIn("NOTE: Returning to sequence selection.", s.output)
        self.assertIsNotNone(scenarios()["main_protein_exit"].exit)       # [0] Exit leaves Synapse

    def test_sequence_selector(self):
        s = scenarios()["selector"]
        self.assertEqual(s.menus[0]["entries"], [("1", "a  (90 bp, linear)"), ("2", "b  (60 bp, circular)"),
                                                 ("F", "Load another file")])
        self.assertIsNone(s.menus[0]["zero"])
        self.assertEqual(s.output.count(style.INVALID_SELECTION), 2)
        self.assertEqual(s.result[0], 2)
        self.assertEqual(scenarios()["selector_file"].result, (None, None))


class MissingFileTests(unittest.TestCase):
    """Master Prompt 11.4: a missing file is reported and the dialog stays usable."""

    MISSING = str(PROJECT_ROOT / "no_such_dir" / "missing.fasta")

    def test_loaders_raise_a_clear_error(self):
        from analyses import blast, crispr, variant_analysis
        from fileio.loader import load_sequences
        cases = [
            (lambda: load_sequences(self.MISSING), FileNotFoundError, "Sequence file not found"),
            (lambda: load_sequences(""), FileNotFoundError, "Sequence file not found"),
            (lambda: crispr.load_fasta(self.MISSING), crispr.CRISPRAnalysisError, "Unable to read FASTA"),
            (lambda: blast.load_fasta_queries(self.MISSING), blast.BlastAnalysisError, "FASTA file not found"),
            (lambda: variant_analysis.parse_vcf_file(self.MISSING), FileNotFoundError, "missing.fasta"),
        ]
        for call, error, text in cases:
            with self.subTest(error=error.__name__, text=text):
                with self.assertRaises(error) as raised:
                    call()
                self.assertIn(text, str(raised.exception))

    def test_menus_report_the_missing_file_and_continue(self):
        from analyses import blast, crispr, variant_analysis
        manager = blast.LocalBlastDatabaseManager(str(PROJECT_ROOT / "no_such_dir"))
        cases = [
            ("variant", ["1", self.MISSING, "0"], variant_analysis_menu.run_variant_analysis_ui,
             "ERROR: Could not read VCF file:"),
            ("crispr batch", ["2", "1", self.MISSING, "0", "0"], lambda: crispr_menu.crispr_menu(DNA, "demo"),
             "ERROR: Unable to read FASTA:"),
            ("blast query", ["1", "2", self.MISSING, "", "0"], lambda: blast_menu.blast_menu(DNA, "demo"),
             "ERROR: FASTA file not found:"),
            ("blast database", ["4", self.MISSING, "db", "1", "0"],
             lambda: blast_menu.blast_database_manager_menu(manager), "ERROR: FASTA file not found:"),
        ]
        for name, answers, call, message in cases:
            with self.subTest(menu=name):
                session = Session(answers).run(call)
                self.assertIn(message, session.output)
                self.assertNotIn("Traceback", session.output)
        self.assertFalse((PROJECT_ROOT / "no_such_dir").exists())   # nothing was created

    def test_main_loader_asks_again_after_a_missing_file(self):
        import main
        with self.assertRaises(_Exhausted):        # it keeps asking; the script ends it
            Session([self.MISSING]).run(main.choose_sequence_file)
        session = Session([self.MISSING, "Q"])
        session.run(main.choose_sequence_file)
        self.assertIn("SEQUENCE FILE ERROR", session.output)
        self.assertIn("ERROR: Sequence file not found:", session.output)


if __name__ == "__main__":
    unittest.main()
