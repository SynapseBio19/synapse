# ============================================================
# SYNAPSE - PIPELINE CORE
# ============================================================
#
# Central pipeline for the five current Synapse modules:
#
# 1. Sequence Overview
# 2. ORF Analysis
# 3. RNA Analysis
# 4. Protein Analysis
# 5. Restriction & Cloning
#
# This file does NOT replace the analysis modules.
# It connects their existing functions and keeps their
# results in one shared context.
#
# ============================================================

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from core.models import AnalysisResult, SequenceRecord, sha256_text
from data.restriction_enzymes import ENZYME_DATABASE_VERSION

logger = logging.getLogger("synapse.pipeline")

# Modules whose results depend on the selected ORF; they are invalidated when
# the ORF selection changes.
ORF_DEPENDENT_MODULES = ("rna", "protein", "crispr", "blast")


# ============================================================
# PIPELINE CONTEXT
# ============================================================

# Every module the pipeline reports on, with its status label. ``completed``
# always has exactly these keys, so ``status()`` has a fixed shape (audit D19).
STATUS_LABELS = {
    "sequence": "Sequence Overview",
    "orf": "ORF Analysis",
    "rna": "RNA Analysis",
    "protein": "Protein Analysis",
    "restriction": "Restriction Analysis",
    "cloning": "Cloning Analysis",
    "primer": "Primer Design",
    "mutation": "Mutation Analysis",
    "variant": "Variant Analysis",
    "crispr": "CRISPR Guide Design",
    "blast": "BLAST Search",
    "literature": "Literature Search",
}
STATUS_MODULES = tuple(STATUS_LABELS)

# The interface between core and whatever displays the pipeline status:
# a callable that receives PipelineContext.status_report(). core imports no
# interface module; the terminal passes ui.report.print_pipeline_status.
StatusRenderer = Callable[[Dict[str, Any]], None]


@dataclass
class PipelineContext:
    """
    Shared state for one selected DNA sequence.

    The context is the memory of the current Synapse session.
    """

    # --------------------------------------------------------
    # SOURCE SEQUENCE
    # --------------------------------------------------------

    sequence: str = ""
    header: str = ""
    sequence_id: Optional[int] = None
    topology: str = "linear"
    record: Optional[SequenceRecord] = None

    # Generic result store: module key -> list of AnalysisResult (newest last).
    results: Dict[str, List[AnalysisResult]] = field(default_factory=dict)

    # --------------------------------------------------------
    # ANALYSIS RESULTS
    # --------------------------------------------------------

    sequence_analysis: Optional[Any] = None

    orf_analysis: Optional[Any] = None
    selected_orf: Optional[Any] = None

    rna_analysis: Optional[Any] = None
    rna_sequence: Optional[str] = None

    protein_analysis: Optional[Any] = None
    protein_sequence: Optional[str] = None

    restriction_analysis: Optional[Any] = None

    cloning_analysis: Optional[Any] = None
    cloning_planner: Optional[Any] = None

    # --------------------------------------------------------
    # SOURCE TRACKING
    # --------------------------------------------------------

    protein_source: Optional[str] = None
    rna_source: Optional[str] = None

    # --------------------------------------------------------
    # COMPLETION STATE
    # --------------------------------------------------------

    completed: dict = field(
        default_factory=lambda: dict.fromkeys(STATUS_MODULES, False)
    )

    # ========================================================
    # LOAD / RESET SEQUENCE
    # ========================================================

    def set_sequence(
        self,
        sequence: str,
        header: str = "Unnamed Sequence",
        sequence_id: Optional[int] = None,
        topology: str = "linear",
    ):
        """Load a new DNA sequence and invalidate old results."""

        self.sequence = (
            sequence
            .replace(" ", "")
            .replace("\n", "")
            .replace("\r", "")
            .upper()
        )

        self.header = header
        self.sequence_id = sequence_id
        self.topology = topology
        self.record = None

        self.clear_results()

    # ========================================================
    # CLEAR ANALYSIS RESULTS
    # ========================================================

    def clear_results(self):
        """Clear every result derived from the current sequence."""

        self.sequence_analysis = None

        self.orf_analysis = None
        self.selected_orf = None

        self.rna_analysis = None
        self.rna_sequence = None

        self.protein_analysis = None
        self.protein_sequence = None

        self.restriction_analysis = None

        self.cloning_analysis = None
        self.cloning_planner = None

        self.protein_source = None
        self.rna_source = None

        self.results = {}

        self.completed = dict.fromkeys(STATUS_MODULES, False)

    # ========================================================
    # BASIC INFORMATION
    # ========================================================

    @property
    def length(self) -> int:
        return len(self.sequence)

    def has_sequence(self) -> bool:
        return bool(self.sequence)

    # ========================================================
    # SEQUENCE
    # ========================================================

    # ========================================================
    # GENERIC RESULT STORE
    # ========================================================

    @property
    def sequence_sha256(self) -> Optional[str]:
        return sha256_text(self.sequence) if self.sequence else None

    def store_result(self, result: AnalysisResult) -> AnalysisResult:
        """Store an AnalysisResult under its module key (newest last)."""
        if result.input_sequence_sha256 is None:
            result.input_sequence_sha256 = self.sequence_sha256
        if not result.input_description:
            result.input_description = self.header
        self.results.setdefault(result.module, []).append(result)
        self.completed[result.module] = True
        logger.info("Stored %s result (%d total)", result.module, len(self.results[result.module]))
        return result

    def store(
        self,
        module: str,
        data: Any,
        parameters: Optional[Dict[str, Any]] = None,
        warnings: Optional[List[str]] = None,
        limitations: Optional[List[str]] = None,
    ) -> AnalysisResult:
        """Wrap ``data`` in an AnalysisResult and store it."""
        return self.store_result(
            AnalysisResult(
                module=str(module),
                data=data,
                parameters=dict(parameters or {}),
                warnings=list(warnings or []),
                limitations=list(limitations or []),
            )
        )

    def latest(self, module: str) -> Optional[AnalysisResult]:
        entries = self.results.get(module) or []
        return entries[-1] if entries else None

    # ========================================================
    # SEQUENCE
    # ========================================================

    def store_sequence_analysis(
        self,
        result: Any
    ):
        self.sequence_analysis = result
        self.completed["sequence"] = True
        self.store("sequence", result)

    # ========================================================
    # ORF
    # ========================================================

    def store_orf_analysis(
        self,
        result: Any
    ):
        self.orf_analysis = result
        self.completed["orf"] = True
        self.store("orf", result, parameters=(result or {}).get("parameters", {}) if isinstance(result, dict) else {})

    def set_selected_orf(
        self,
        orf: Optional[Any]
    ):
        """
        Store the selected ORF.

        Changing the ORF invalidates RNA and protein because
        both can depend on the selected coding sequence.
        """

        self.selected_orf = orf

        self.rna_analysis = None
        self.rna_sequence = None

        self.protein_analysis = None
        self.protein_sequence = None

        self.rna_source = None
        self.protein_source = None

        self.completed["rna"] = False
        self.completed["protein"] = False

        for module in ORF_DEPENDENT_MODULES:
            if self.results.pop(module, None):
                logger.info("ORF selection changed: cleared stored %s results", module)

    # ========================================================
    # RNA
    # ========================================================

    def store_rna_analysis(
        self,
        result: Any,
        rna_sequence: Optional[str] = None
    ):
        self.rna_analysis = result
        self.rna_sequence = rna_sequence
        self.completed["rna"] = True
        source = (result or {}).get("pipeline_input", {}).get("source") if isinstance(result, dict) else None
        self.store("rna", result, parameters={"source": source or self.rna_source})

    # ========================================================
    # PROTEIN
    # ========================================================

    def store_protein_analysis(
        self,
        result: Any,
        protein_sequence: Optional[str] = None
    ):
        self.protein_analysis = result
        self.protein_sequence = protein_sequence
        self.completed["protein"] = True
        self.store("protein", result, parameters={"source": self.protein_source})

    # ========================================================
    # RESTRICTION
    # ========================================================

    def store_restriction_analysis(
        self,
        result: Any
    ):
        self.restriction_analysis = result
        self.completed["restriction"] = True
        parameters: Dict[str, Any] = {"enzyme_database_version": ENZYME_DATABASE_VERSION}
        if isinstance(result, dict):
            parameters["enzymes"] = list(result.get("enzymes") or [])
            parameters["topology"] = result.get("topology")
        self.store("restriction", result, parameters=parameters)

    # ========================================================
    # CLONING
    # ========================================================

    def store_cloning_analysis(
        self,
        result: Any
    ):
        self.cloning_analysis = result
        self.completed["cloning"] = True
        self.store("cloning", result,
                   parameters={"enzyme_database_version": ENZYME_DATABASE_VERSION})

    def store_cloning_planner(
        self,
        result: Any
    ):
        self.cloning_planner = result
        self.completed["cloning"] = True
        self.store("cloning", result, parameters={"mode": "automatic planner",
                                                "enzyme_database_version": ENZYME_DATABASE_VERSION})

    # ========================================================
    # PRIMER (called by analyses.primer UI)
    # ========================================================

    def store_primer_analysis(
        self,
        result: Any
    ):
        # A design around a restriction site (workflow 10.2) wraps the
        # primer report; its site and search regions are recorded too.
        site_design = hasattr(result, "enzyme") and hasattr(result, "report")
        source = result.report if site_design else result
        parameters = {}
        for name in ("tm_method", "target_start", "target_end"):
            if hasattr(source, name):
                parameters[name] = getattr(source, name)
        conditions = getattr(source, "thermodynamic_conditions", None)
        if conditions is not None:
            try:
                from dataclasses import asdict as _asdict
                parameters["thermodynamic_conditions"] = _asdict(conditions)
            except TypeError:
                parameters["thermodynamic_conditions"] = conditions
        warnings = list(getattr(source, "warnings", []) or [])
        if site_design:
            for name in ("enzyme", "site_start", "site_end", "flank"):
                parameters[name] = getattr(result, name)
            parameters["forward_region"] = list(result.forward_region)
            parameters["reverse_region"] = list(result.reverse_region)
            warnings += [w for w in result.warnings if w not in warnings]
        self.store("primer", result, parameters=parameters, warnings=warnings)

    # ========================================================
    # STATUS
    # ========================================================

    def status(self):
        return self.completed.copy()

    # ========================================================
    # SUMMARY
    # ========================================================

    def summary(self):

        return {
            "sequence": self.sequence,
            "header": self.header,
            "sequence_id": self.sequence_id,
            "topology": self.topology,

            "sequence_analysis":
                self.sequence_analysis,

            "orf_analysis":
                self.orf_analysis,

            "selected_orf":
                self.selected_orf,

            "rna_analysis":
                self.rna_analysis,

            "rna_sequence":
                self.rna_sequence,

            "protein_analysis":
                self.protein_analysis,

            "protein_sequence":
                self.protein_sequence,

            "restriction_analysis":
                self.restriction_analysis,

            "cloning_analysis":
                self.cloning_analysis,

            "cloning_planner":
                self.cloning_planner,

            "rna_source":
                self.rna_source,

            "protein_source":
                self.protein_source,

            "completed":
                self.completed.copy(),
        }

    # ========================================================
    # STATUS REPORT
    # ========================================================

    def status_report(self) -> Dict[str, Any]:
        """The pipeline status as plain data, for whichever interface shows it.

        core decides what the status is; how it is drawn belongs to the
        interface (``ui.report.print_pipeline_status`` for the terminal).
        """
        return {
            "header": self.header,
            "length": self.length,
            "topology": self.topology,
            "sequence_sha256": self.sequence_sha256 if self.sequence else None,
            "source_path": (self.record.source_path
                            if self.record is not None and self.record.source_path else None),
            # None: no ORF selected; "-": an ORF without a rank.
            "selected_orf_rank": (None if self.selected_orf is None
                                  else self.selected_orf.get("rank", "-")),
            "rna_source": self.rna_source or None,
            "protein_source": self.protein_source or None,
            "analyses": [
                {"key": key, "label": label, "stored": len(self.results.get(key, []))}
                for key, label in STATUS_LABELS.items()
            ],
        }

    def print_status(self, render: StatusRenderer) -> None:
        """Hand the status to an interface's renderer; core itself draws nothing."""
        render(self.status_report())


# ============================================================
# SYNAPSE PIPELINE
# ============================================================

class SynapsePipeline:
    """
    High-level controller for the current Synapse pipeline.

    Existing analysis modules remain responsible for the actual
    calculations. This class only connects them and stores results.
    """

    def __init__(
        self,
        sequence: str = "",
        header: str = "",
        sequence_id: Optional[int] = None,
        topology: str = "linear",
    ):

        self.context = PipelineContext()

        if sequence:

            self.context.set_sequence(
                sequence=sequence,
                header=header,
                sequence_id=sequence_id,
                topology=topology,
            )

    # ========================================================
    # RESULT STORE (delegates to the context)
    # ========================================================

    def attach_record(self, record: SequenceRecord) -> None:
        """Attach provenance (file, hash, topology, sequence type) for the loaded sequence.

        An RNA record matches a loaded sequence in the DNA alphabet (U read as T).
        """
        expected = record.sequence
        if record.alphabet == "RNA":
            expected = expected.replace("U", "T")
        if self.context.sequence not in (record.sequence, expected):
            raise ValueError("SequenceRecord does not match the loaded sequence.")
        self.context.record = record
        self.context.topology = record.topology

    def store(self, module: str, data: Any, **kwargs) -> AnalysisResult:
        return self.context.store(module, data, **kwargs)

    def store_result(self, result: AnalysisResult) -> AnalysisResult:
        return self.context.store_result(result)

    @property
    def results(self) -> Dict[str, List[AnalysisResult]]:
        return self.context.results

    def latest(self, module: str) -> Optional[AnalysisResult]:
        return self.context.latest(module)

    # ========================================================
    # LOAD SEQUENCE
    # ========================================================

    def load_sequence(
        self,
        sequence: str,
        header: str = "Unnamed Sequence",
        sequence_id: Optional[int] = None,
        topology: str = "linear",
    ):

        self.context.set_sequence(
            sequence=sequence,
            header=header,
            sequence_id=sequence_id,
            topology=topology,
        )

    # ========================================================
    # SEQUENCE ANALYSIS
    # ========================================================

    def run_sequence(
        self,
        **kwargs
    ):
        """
        Run the existing analyses.sequence.analyze_sequence().
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No DNA sequence loaded."
            )

        from analyses.sequence import (
            analyze_sequence
        )

        result = analyze_sequence(
            self.context.sequence,
            **kwargs
        )

        self.context.store_sequence_analysis(
            result
        )

        return result

    # ========================================================
    # ORF ANALYSIS
    # ========================================================

    def run_orf(
        self,
        include_partial=True,
        min_length=None,
        enhanced=True,
        top_n=5,
        start_codon_policy=None,
        nested=None,
        ranking=None,
        genetic_code=1,
    ):
        """
        Run the existing ORF analysis.

        Uses analyze_orfs_enhanced() when available.
        Falls back to analyze_orfs() if necessary.

        start_codon_policy, nested, ranking and genetic_code are the ORF
        menu's search settings; None means the module default, so a
        programmatic run can reproduce an interactive one (audit D19).
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No DNA sequence loaded."
            )

        from analyses import orf as orf_module

        if min_length is None:

            min_length = getattr(
                orf_module,
                "DEFAULT_MIN_ORF_LENGTH",
                10
            )

        search_settings = {
            "start_codon_policy": start_codon_policy,
            "nested": orf_module.DEFAULT_REPORT_NESTED if nested is None else nested,
            "ranking": ranking or orf_module.DEFAULT_RANKING,
            "genetic_code": genetic_code,
            "topology": self.context.topology,
        }

        if (
            enhanced
            and hasattr(
                orf_module,
                "analyze_orfs_enhanced"
            )
        ):

            result = (
                orf_module.analyze_orfs_enhanced(
                    self.context.sequence,
                    include_partial=include_partial,
                    min_length=min_length,
                    top_n=top_n,
                    **search_settings,
                )
            )

        else:

            result = (
                orf_module.analyze_orfs(
                    self.context.sequence,
                    include_partial=include_partial,
                    min_length=min_length,
                    **search_settings,
                )
            )

        self.context.store_orf_analysis(
            result
        )

        return result

    # ========================================================
    # SELECT ORF
    # ========================================================

    def set_selected_orf(
        self,
        orf: Optional[Any]
    ):
        """
        Store an already selected ORF.

        This is the preferred method when the UI has already
        handled ORF selection.
        """

        self.context.set_selected_orf(
            orf
        )

        return orf

    def select_orf(
        self,
        interactive=True
    ):
        """
        Select an ORF from the current ORF analysis.

        interactive=True uses the existing select_orf() UI.
        """

        if self.context.orf_analysis is None:

            self.run_orf()

        orf_result = (
            self.context.orf_analysis
        )

        if isinstance(
            orf_result,
            dict
        ):

            orfs = orf_result.get(
                "orfs",
                []
            )

        else:

            orfs = orf_result

        if not orfs:

            self.context.set_selected_orf(
                None
            )

            return None

        if interactive:

            from analyses.orf import (
                select_orf
            )

            selected = select_orf(
                orfs
            )

        else:

            selected = orfs[0]

        self.context.set_selected_orf(
            selected
        )

        return selected

    # ========================================================
    # RNA SOURCE
    # ========================================================

    def get_rna_source_dna(self):
        """
        Return DNA for RNA analysis.

        Priority:
            1. selected ORF
            2. complete selected DNA
        """

        if self.context.selected_orf is not None:

            sequence = (
                self.context.selected_orf.get(
                    "sequence"
                )
            )

            if sequence:

                return (
                    sequence,
                    "Selected ORF"
                )

        return (
            self.context.sequence,
            "Selected DNA"
        )

    # ========================================================
    # RNA ANALYSIS
    # ========================================================

    def run_rna(
        self,
        use_selected_orf=True,
        **kwargs
    ):
        """
        Run the existing analyses.rna.analyze_rna().

        If an ORF is selected, its coding DNA is used.
        Otherwise the selected DNA sequence is used.
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No DNA sequence loaded."
            )

        from analyses.rna import (
            analyze_rna,
            transcribe_rna,
        )

        if use_selected_orf:

            dna_source, source_name = (
                self.get_rna_source_dna()
            )

        else:

            dna_source = (
                self.context.sequence
            )

            source_name = (
                "Selected DNA"
            )

        rna = transcribe_rna(
            dna_source
        )

        # Audit P1-9: a selected ORF carries the genetic code it was found
        # with, so the RNA codon view uses it unless the caller says otherwise.
        if use_selected_orf and self.context.selected_orf is not None:
            from analyses.rna import orf_genetic_code

            kwargs.setdefault(
                "genetic_code",
                orf_genetic_code(self.context.selected_orf),
            )

        result = analyze_rna(
            dna_source,
            **kwargs
        )

        self.context.rna_source = (
            source_name
        )

        self.context.store_rna_analysis(
            result,
            rna_sequence=rna
        )

        return result

    # ========================================================
    # PROTEIN SOURCE
    # ========================================================

    def get_protein_source_rna(self):
        """
        Return RNA for protein analysis.

        If a selected ORF exists, its coding sequence is used.
        Otherwise the current DNA sequence is transcribed.
        """

        if self.context.selected_orf is not None:

            if self.context.rna_sequence:

                return (
                    self.context.rna_sequence,
                    "Selected ORF"
                )

            from analyses.rna import (
                transcribe_rna
            )

            dna = (
                self.context.selected_orf.get(
                    "sequence",
                    ""
                )
            )

            return (
                transcribe_rna(dna),
                "Selected ORF"
            )

        if self.context.rna_sequence:

            return (
                self.context.rna_sequence,
                self.context.rna_source
                or "Selected DNA"
            )

        from analyses.rna import (
            transcribe_rna
        )

        return (
            transcribe_rna(
                self.context.sequence
            ),
            "Selected DNA"
        )

    # ========================================================
    # PROTEIN ANALYSIS
    # ========================================================

    def run_protein(
        self,
        use_selected_orf=True
    ):
        """
        Run the existing complete protein analysis.

        The current protein module accepts:
            analyze_complete_protein(
                protein="",
                rna=<RNA>
            )

        Therefore this pipeline supplies RNA derived from the
        selected ORF when one exists.
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No DNA sequence loaded."
            )

        from analyses.protein import (
            analyze_complete_protein
        )

        if (
            use_selected_orf
            and self.context.selected_orf is not None
        ):

            rna, source_name = (
                self.get_protein_source_rna()
            )
            # A 5'-partial ORF has no start codon in the sequence.
            initiator_start = (source_name == "Selected ORF"
                               and self.context.selected_orf.get("status") not in (
                                   "Partial - Missing Start", "Partial - Missing Start and Stop"))

        else:

            from analyses.rna import (
                transcribe_rna
            )

            rna = transcribe_rna(
                self.context.sequence
            )

            source_name = (
                "Selected DNA"
            )
            initiator_start = False

        result = analyze_complete_protein(
            protein="",
            rna=rna,
            initiator_start=initiator_start,
            genetic_code=(
                self.context.selected_orf.get("genetic_code", 1)
                if initiator_start and isinstance(self.context.selected_orf, dict)
                else 1
            ),
        )

        protein = result.get(
            "protein",
            ""
        )

        self.context.protein_source = (
            source_name
        )

        self.context.store_protein_analysis(
            result,
            protein_sequence=protein
        )

        return result

    # ========================================================
    # RESTRICTION ANALYSIS
    # ========================================================

    def run_restriction(
        self,
        enzymes,
        topology=None
    ):
        """
        Run the existing restriction analyzer.

        Parameters:
            enzymes:
                list of enzyme names

            topology:
                'linear' or 'circular'
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No DNA sequence loaded."
            )

        from analyses.restriction import (
            analyze_restriction
        )

        if topology is None:

            topology = (
                self.context.topology
            )

        result = analyze_restriction(
            self.context.sequence,
            enzymes,
            topology
        )

        self.context.store_restriction_analysis(
            result
        )

        return result

    # ========================================================
    # CLONING ANALYSIS
    # ========================================================

    def run_cloning(
        self,
        insert_sequence,
        vector_left_enzyme,
        vector_right_enzyme,
        insert_left_enzyme,
        insert_right_enzyme,
    ):
        """
        Run the existing two-enzyme cloning analyzer.

        The current selected DNA is always the vector.
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No vector DNA sequence loaded."
            )

        from analyses.cloning import (
            analyze_enzyme_pair
        )

        result = analyze_enzyme_pair(
            self.context.sequence,
            insert_sequence,
            vector_left_enzyme,
            vector_right_enzyme,
            insert_left_enzyme,
            insert_right_enzyme,
            vector_topology=self.context.topology,
        )

        self.context.store_cloning_analysis(
            result
        )

        return result

    # ========================================================
    # AUTOMATIC CLONING PLANNER
    # ========================================================

    def run_cloning_planner(
        self,
        insert_sequence,
        enzyme_names=None,
        top_n=5
    ):
        """
        Run the existing automatic cloning planner.

        The current selected DNA is used as the vector.
        """

        if not self.context.has_sequence():

            raise ValueError(
                "No vector DNA sequence loaded."
            )

        from analyses.cloning import (
            automatic_cloning_planner
        )

        result = automatic_cloning_planner(
            self.context.sequence,
            insert_sequence,
            enzyme_names,
            top_n,
            vector_topology=self.context.topology,
        )

        self.context.store_cloning_planner(
            result
        )

        return result

    # ========================================================
    # RESTRICTION + CLONING
    # ========================================================

    def run_restriction_and_cloning(
        self,
        restriction_enzymes,
        insert_sequence,
        vector_left_enzyme,
        vector_right_enzyme,
        insert_left_enzyme,
        insert_right_enzyme,
        topology=None,
    ):
        """
        Convenience method for a complete restriction/cloning
        workflow using the same selected vector DNA.
        """

        restriction_result = (
            self.run_restriction(
                restriction_enzymes,
                topology=topology
            )
        )

        cloning_result = (
            self.run_cloning(
                insert_sequence,
                vector_left_enzyme,
                vector_right_enzyme,
                insert_left_enzyme,
                insert_right_enzyme,
            )
        )

        return {
            "restriction":
                restriction_result,

            "cloning":
                cloning_result,
        }

    # ========================================================
    # COMPLETE CORE PIPELINE
    # ========================================================

    def run_core_pipeline(
        self,
        run_sequence=True,
        run_orf=True,
        run_rna=True,
        run_protein=True,
        restriction_enzymes=None,
        restriction_topology=None,
        include_partial=True,
        min_orf_length=None,
        top_n=5,
    ):
        """
        Run the sequence -> ORF -> RNA -> protein chain.

        Restriction is optional because enzyme selection is a
        user/design decision.

        Cloning is intentionally not run here because it requires
        an insert and cloning enzyme choices.
        """

        results = {}

        # ----------------------------------------------------
        # SEQUENCE
        # ----------------------------------------------------

        if run_sequence:

            results["sequence"] = (
                self.run_sequence()
            )

        # ----------------------------------------------------
        # ORF
        # ----------------------------------------------------

        if run_orf:

            results["orf"] = self.run_orf(
                include_partial=include_partial,
                min_length=min_orf_length,
                enhanced=True,
                top_n=top_n,
            )

        # ----------------------------------------------------
        # RNA
        # ----------------------------------------------------

        if run_rna:

            results["rna"] = self.run_rna(
                use_selected_orf=True
            )

        # ----------------------------------------------------
        # PROTEIN
        # ----------------------------------------------------

        if run_protein:

            results["protein"] = (
                self.run_protein(
                    use_selected_orf=True
                )
            )

        # ----------------------------------------------------
        # RESTRICTION
        # ----------------------------------------------------

        if restriction_enzymes:

            results["restriction"] = (
                self.run_restriction(
                    restriction_enzymes,
                    topology=restriction_topology
                )
            )

        return results

    # ========================================================
    # PIPELINE STATE
    # ========================================================

    def status(self):

        return self.context.status()

    def status_report(self) -> Dict[str, Any]:

        return self.context.status_report()

    def print_status(self, render: StatusRenderer) -> None:

        self.context.print_status(render)

    def get_results(self):

        return self.context.summary()


# ============================================================
# FACTORY
# ============================================================

def create_pipeline(
    sequence: str,
    header: str = "Unnamed Sequence",
    sequence_id: Optional[int] = None,
    topology: str = "linear",
):

    return SynapsePipeline(
        sequence=sequence,
        header=header,
        sequence_id=sequence_id,
        topology=topology,
    )


# ============================================================
# TEST
# ============================================================

if __name__ == "__main__":
    # A smoke test of core alone: it prints the status as data, so it runs
    # without any interface module.
    import json

    test_sequence = (
        "ATGCGTACGAATTCGGGATCCATGCGAAGCTTACGTCGAGCTCGGTACCT"
        "CTAGAGGATCCGATATCGCGTACG"
    )

    pipeline = create_pipeline(
        sequence=test_sequence,
        header="TEST_DNA",
        sequence_id=1,
        topology="linear",
    )

    print(json.dumps(pipeline.status_report(), indent=2))