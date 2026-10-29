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

from dataclasses import dataclass, field
from typing import Any, Optional


# ============================================================
# PIPELINE CONTEXT
# ============================================================

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
        default_factory=lambda: {
            "sequence": False,
            "orf": False,
            "rna": False,
            "protein": False,
            "restriction": False,
            "cloning": False,
        }
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

        self.completed = {
            "sequence": False,
            "orf": False,
            "rna": False,
            "protein": False,
            "restriction": False,
            "cloning": False,
        }

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

    def store_sequence_analysis(
        self,
        result: Any
    ):
        self.sequence_analysis = result
        self.completed["sequence"] = True

    # ========================================================
    # ORF
    # ========================================================

    def store_orf_analysis(
        self,
        result: Any
    ):
        self.orf_analysis = result
        self.completed["orf"] = True

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

    # ========================================================
    # RESTRICTION
    # ========================================================

    def store_restriction_analysis(
        self,
        result: Any
    ):
        self.restriction_analysis = result
        self.completed["restriction"] = True

    # ========================================================
    # CLONING
    # ========================================================

    def store_cloning_analysis(
        self,
        result: Any
    ):
        self.cloning_analysis = result
        self.completed["cloning"] = True

    def store_cloning_planner(
        self,
        result: Any
    ):
        self.cloning_planner = result
        self.completed["cloning"] = True

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

    def print_status(self):

        labels = {
            "sequence": "Sequence Overview",
            "orf": "ORF Analysis",
            "rna": "RNA Analysis",
            "protein": "Protein Analysis",
            "restriction": "Restriction Analysis",
            "cloning": "Cloning Analysis",
        }

        print()
        print("=" * 70)
        print("SYNAPSE PIPELINE STATUS")
        print("=" * 70)

        print()
        print("Current Sequence")
        print("-" * 45)

        print(
            f"Header:       "
            f"{self.header or '-'}"
        )

        print(
            f"Length:       "
            f"{self.length} bp"
        )

        print(
            f"Topology:     "
            f"{self.topology}"
        )

        if self.selected_orf is not None:

            print(
                f"Selected ORF: "
                f"#{self.selected_orf.get('rank', '-')}"
            )

        if self.rna_source:

            print(
                f"RNA Source:   "
                f"{self.rna_source}"
            )

        if self.protein_source:

            print(
                f"Protein Source:"
                f" {self.protein_source}"
            )

        print()
        print("Analysis Status")
        print("-" * 45)

        for key, label in labels.items():

            state = (
                "COMPLETED"
                if self.completed[key]
                else "NOT RUN"
            )

            print(
                f"{label:<28} "
                f"{state}"
            )

        print()
        print("=" * 70)


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
    ):
        """
        Run the existing ORF analysis.

        Uses analyze_orfs_enhanced() when available.
        Falls back to analyze_orfs() if necessary.
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
                )
            )

        else:

            result = (
                orf_module.analyze_orfs(
                    self.context.sequence,
                    include_partial=include_partial,
                    min_length=min_length,
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
            initiator_start = source_name == "Selected ORF"

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
            top_n
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

    def print_status(self):

        self.context.print_status()

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

    print()
    print("=" * 70)
    print("SYNAPSE PIPELINE TEST")
    print("=" * 70)

    print(
        f"\nSequence: "
        f"{pipeline.context.header}"
    )

    print(
        f"Length:   "
        f"{pipeline.context.length} bp"
    )

    print(
        "\nPipeline initialized successfully."
    )

    pipeline.print_status()