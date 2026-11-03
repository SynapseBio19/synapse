# Feature verification matrix

State as of the Phase 6 initial audit, 2026-10-30, commit `0dad3bb`,
updated through the end of Phase 6 (2026-11-02). A finding marked
corrected names the date; the evidence is in `docs/CHANGES.md` and
`docs/PHASE6_PROGRESS.md`.
Menu numbers are those of the code at the end of Stage F (2026-11-03);
the line and definition counts in the section headings are those of the
initial audit (non-blank lines at `0dad3bb`).
Companion to `docs/PHASE6_INITIAL_AUDIT.md`, which holds the evidence for
every defect referenced here as `D…` (confirmed) or `S…` (suspected).

The audit itself modified no application source file; later rows record
the corrections made since.

## Status vocabulary

| Status | Meaning |
|---|---|
| **Implemented** | The code exists and is reachable. Says nothing about correctness. |
| **Execution-tested** | Actually run in this audit or by an automated test, with its output inspected. |
| **Integration-tested** | Exercised together with at least one dependent module, including the hand-over of data. |
| **Independently validated** | Compared against an external reference implementation or a published worked example. Only used where that comparison was really performed. |
| **Partially implemented** | Some advertised part is missing; the missing part is named. |
| **Not implemented** | Advertised somewhere (menu, docs, config) but absent. |
| **Blocked by external dependency** | Cannot be completed or verified here; the blocker is named. |
| **Not applicable** | The column does not apply to this item. |

"Reachable" in the execution column means the menu entry dispatches to real
code that was read, but the interactive path itself was not driven by
simulated keystrokes (see audit §11.5).

## Summary

| | Count |
|---|---|
| Analysis modules | 12 |
| Core modules | 9 (`core/alphabets.py` added in Phase 6) |
| Shared UI module | 1 (`ui/style.py`, added in Stage F) |
| Support modules (`fileio`, `ui`, `data`) | 8 |
| Functions and classes inventoried | ~1,100 |
| Main-menu entries | 15 |
| Sub-menu entries inventoried | 78 |
| Automated tests | 356, all passing (2026-11-03). The count at the end of Stage 4 was 283: (109 at the audit + 9 D12 + 7 D4 + 7 D3 + 7 D5 + 7 D1 + 5 D2 + 5 D14-F2 + 11 sequence-type classifier + 6 D14-F3/F11 + 7 D14-F1 + 9 D13/D14-F4/F5/F7 + 5 D14-F6 + 3 D14-F8–F10 + 3 D6 + 3 D9 + 3 D8 + 4 D10/D11 + 2 D19 + 1 D15 + 2 S15/S16 + 2 S17 + 2 S19 + 9 S23/S3/S4/S5 + 1 D16 + 7 S8/S10/S21/S22 + 4 S1/S2/S9 + 5 P1-9 + 3 S18 + 8 S14 + 9 P1-7 + 2 S13 + 12 Stage F + 4 S24) |
| Module self-tests | 12 suites, all passing (literature 64, mutation 44, BLAST 33, CRISPR 28, variant 27, others aggregated) |
| Golden-master scenarios | 12 |
| Confirmed defects | 19 (D1–D19), all corrected by the end of Phase 6 (`docs/PHASE6_FINAL_REPORT.md`, section A; D16 finished with the removal of the superseded primer generations). The first 16 corrections are dated 2026-10-31 in the rows below. |
| Suspected / unverified concerns | 14 at the audit (S1–S14); S15–S24 recorded since |
| External validations performed in this audit | 8 (audit §3, §6) |
| Features advertised but not implemented | 1 (`[11] AI Experiment Planner`, labeled "not implemented" since D15) |

---

## 1. Main menu — `main.py`, `ui/menu.py`

| # | Feature | Function | Implementation | Execution | Integration | Validation | Notes |
|---|---|---|---|---|---|---|---|
| — | Load FASTA/GenBank file | `main.choose_sequence_file`, `fileio.loader.load_sequences` | Implemented | Execution-tested | Integration-tested | Not applicable | Rejects empty and header-less files; content-sniffs, ignores extension. **D14 corrected 2026-10-31**: each record is classified by `core.alphabets`; a record whose type the letters cannot decide is not guessed, the user declares it. |
| — | Sequence selection, `[F]` load another file | `ui.menu.select_gene` | Implemented | Execution-tested | Integration-tested | Not applicable | Shows length and topology. `examples/dna.fasta` = 14 records. |
| — | Sequence-type gate | `main.resolve_sequence_type`, `core.validation.validate_dna`, `get_dna_errors` | Implemented | Execution-tested | Integration-tested | Not applicable | Type from `core.alphabets`; undetermined records are declared by the user; RNA analyzed in the DNA alphabet with a note (**D13**, **D14-F4/F5** corrected 2026-10-31). |
| 1 | Sequence Overview | `analyses.sequence.*` (15 functions called from `main.py`) | Implemented | Execution-tested | Integration-tested | Not applicable | Self-test passes. Stored in the pipeline. |
| 2 | ORF Analysis | `analyses.orf.run_orf_analysis_ui` | Implemented | Execution-tested | Integration-tested | Independently validated (genetic codes, coordinates) | See §3. |
| 3 | RNA Analysis | `analyses.rna.run_rna_analysis_ui` | Implemented | Execution-tested | Integration-tested | Not applicable | Selected-ORF vs whole-sequence rule confirmed intact. |
| 4 | Protein Analysis | `analyses.protein.run_protein_analysis_from_orf` | Implemented | Execution-tested | Integration-tested | Independently validated (ProtParam) | Requires a selected ORF by design; no silent whole-sequence fallback. |
| — | Protein record menu | `main.run_protein_record`, `analyses.protein.run_protein_analysis_from_record` | Implemented | Execution-tested | Integration-tested | Not applicable | Protein FASTA/GenPept records: protein analysis, BLAST (protein query), report (D14, 2026-10-31). X/B/Z/J/U/O keep their positions (**D14-F6** corrected 2026-10-31). |
| 5 | Restriction & Cloning | `analyses.restriction.run_restriction_cloning_ui` | Implemented | Execution-tested | Integration-tested | Independently validated (digests, circular map) | D3/D4/D5 corrected 2026-10-31. Open: **S17**, **S18**. |
| 6 | Primer Design | `analyses.primer.run_primer_analysis_ui` | Implemented | Execution-tested | Integration-tested | Independently validated (primer3-py) | **D16 corrected 2026-10-31**: one definition per name; superseded generations removed after a reachability review. Tm agrees with primer3-py within 0.0032 °C over 2,916 comparisons. |
| 7 | Mutation & Variant Analysis | `analyses.mutation.run_mutation_variant_ui` | Implemented | Execution-tested | Integration-tested | Partially (algorithm checked, no bcftools) | **D12** in the variant branch — corrected 2026-10-31 (Phase 6 Tier 1). |
| 8 | CRISPR Guide Design | `analyses.crispr.crispr_menu` | Partially implemented | Execution-tested | Integration-tested | Not applicable | **D8, D9, D10 corrected 2026-10-31** (cut site, on-target exclusion, configuration fields). Partial because CFD scoring and genome-wide off-target search are absent by design. |
| 9 | BLAST Search | `analyses.blast.blast_menu` | Implemented | Execution-tested | Integration-tested | Independently validated (real BLAST+ run) | **D6, D7 corrected 2026-10-31**: programmatic local BLAST by database path; the stale artifact is labeled legacy. |
| 10 | Literature Search | `analyses.literature.main` | Implemented | Execution-tested (offline) | Integration-tested | Blocked by external dependency (no live call made) | **S11**: a live check needs the owner's approval for external requests. |
| 11 | **AI Experiment Planner** | — | **Not implemented** | Not applicable | Not applicable | Not applicable | **D15 corrected 2026-10-31**: the entry says "not implemented" instead of "coming soon". The planner itself is P3-1, deferred. |
| P | Pipeline Status & Stored Results | `core.pipeline.PipelineContext.status_report`, `ui.report.print_pipeline_status` | Implemented | Execution-tested | Integration-tested | Not applicable | **D19 corrected 2026-10-31**: `completed` has a fixed shape. |
| R | Export Analysis Report (JSON/Markdown) | `core.report.export_analysis_report` | Implemented | Execution-tested | Integration-tested | Not applicable | Rejects suffixes other than `.json`/`.md`. **S8 corrected 2026-10-31**: a failing `to_dict()` is logged. Since 2026-11-02 a result's own warnings reach the summary. |
| S | Select Another Sequence | `main.py` loop | Implemented | Execution-tested | Not applicable | Not applicable | |
| 0 | Exit | `main.py` | Implemented | Execution-tested | Not applicable | Not applicable | **D17 corrected 2026-11-01**: `[0]` means Back, Cancel or Exit (main menu only), enforced by `tests/test_terminal_ui.py`. |

---

## 2. Sequence analysis — `analyses/sequence.py` (1,768 lines, 35 functions)

| Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|
| Nucleotide counts, percentages | `count_nucleotides`, `calculate_nucleotide_percentages` | Implemented | Execution-tested | Not applicable | |
| Ambiguous-base counting | `count_ambiguous_bases` | Implemented | Execution-tested | Not applicable | |
| GC / AT content, content ranges | `calculate_gc_content`, `calculate_at_content`, `calculate_dna_content_ranges` | Implemented | Execution-tested | Not applicable | |
| GC / AT skew | `calculate_gc_skew`, `calculate_at_skew` | Implemented | Execution-tested | Not applicable | |
| GC and complexity profiles | `calculate_gc_content_profile`, `calculate_complexity_profile`, `calculate_profile_statistics` | Implemented | Execution-tested | Not applicable | Window/step parameterized via `core.api`. |
| Sequence complexity | `calculate_sequence_complexity` | Implemented | Execution-tested | Not applicable | Normalized Shannon entropy of A/C/G/T composition; defined in `METHODS.md` with its limit (it does not detect repeats). P2-5: the name is kept, the meaning is now stated (2026-11-02). |
| Homopolymer detection | `find_homopolymers` | Implemented | Execution-tested | Not applicable | |
| CpG analysis and islands | `analyze_cpg`, `find_cpg_islands` | Implemented | Execution-tested | Not applicable | P2-6 sliding-count optimization deferred. |
| Molecular weight (DNA) | `calculate_dna_molecular_weight`, `calculate_average_dna_molecular_weight` | Implemented | Execution-tested | Not applicable | |
| Quality assessment | `analyze_sequence_quality` | Implemented | Execution-tested | Not applicable | |
| Complement / reverse complement | `get_complement`, `get_reverse_complement` | Implemented | Execution-tested | Independently validated | Agrees with Biopython on IUPAC symbols. |
| Interpretation text, summary | `create_sequence_interpretation`, `create_sequence_summary` | Implemented | Execution-tested | Not applicable | |
| Unified analysis | `analyze_sequence` | Implemented | Execution-tested | Not applicable | Entry point for `core.api.sequence_overview`. |
| Self-test | `run_sequence_self_test` | Implemented | Execution-tested | Not applicable | All pass. |

Test coverage: **2 test references** — the thinnest coverage of any module
relative to its use (every session runs it).

---

## 3. ORF analysis — `analyses/orf.py` (2,587 lines, 46 functions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| 1 | ORF Discovery | `find_orfs`, `analyze_orfs`, `analyze_orfs_enhanced` | Implemented | Execution-tested | Independently validated (codes, coordinates) | Both strands, 3 frames each. |
| 2 | ORF Candidates | `print_orf_candidates`, `print_orf_candidate_details` | Implemented | Reachable | Not applicable | |
| 3 | Compare ORFs | `compare_orf_candidates`, `print_orf_comparison_table` | Implemented | Execution-tested | Not applicable | Rows carry `heuristic_score` and `heuristic_classification` (**D18 corrected 2026-11-02**). |
| 4 | ORF Relationships | `analyze_orf_relationships` | Implemented | Execution-tested | Not applicable | Overlap/nesting by coordinates only, correctly disclaimed. |
| 5 | Select ORF | `select_orf`, `analyze_selected_orf` | Implemented | Integration-tested | Not applicable | Hands over to RNA, protein, CRISPR, BLAST. |
| 6.1 | Minimum ORF length | `_orf_settings_menu` | Implemented | Execution-tested | Not applicable | Default 75 nt incl. stop; any value ≥ 3. |
| 6.2 | Start-codon policy | `start_codons_for`, `resolve_start_codons` | Implemented | Execution-tested | Independently validated | `atg_only`, `atg_and_alternative`, `table_starts`. |
| 6.3 | Nested ORFs | `_keep_longest_per_stop` | Implemented | Execution-tested | Not applicable | **S2 not a defect** (verified by execution 2026-10-31): nested minus-strand ORFs share frame and 3' end and collapse correctly. |
| 6.4 | Ranking | `_orf_sort_key`, `rank_orf_candidates` | Implemented | Execution-tested | Not applicable | `length` (default) or `heuristic`. |
| 6.5 | Partial ORFs | `find_orfs(include_partial=)` | Partially implemented | Execution-tested | Not applicable | 3′-partial (missing stop) only. **5′-partial not implemented** (METHODS states this). |
| 6.6 | Genetic code | `codon_table_for` | Implemented | Execution-tested | Independently validated | NCBI 1, 2, 4, 11 verified against `Bio.Data.CodonTable`. |
| 6.7/6.8 | Restore Synapse / legacy settings | `LEGACY_ORF_SETTINGS` | Implemented | Execution-tested | Not applicable | Golden-master proves legacy equivalence. |
| — | Start-codon classification | `classify_start_codon` | Implemented | Execution-tested | Not applicable | |
| — | Ambiguity analysis | `analyze_orf_ambiguity` | Implemented | Execution-tested | Not applicable | Deterministic vs uncertain codons separated. |
| — | Codon usage, diversity, GC3 | `analyze_codon_usage`, `calculate_codon_diversity` | Implemented | Execution-tested | Not applicable | |
| — | Kozak context | `analyze_kozak_context` | Implemented | Execution-tested | Not applicable | **S1 not a defect** (verified 2026-10-31): the reverse-complement frame is the biologically correct one. Simplified 3-position heuristic. |
| — | Shine–Dalgarno screen | `analyze_shine_dalgarno` | Implemented | Execution-tested | Not applicable | Fixed motif list, no energy model. |
| — | CDS quality, heuristic score | `calculate_cds_quality`, `calculate_orf_heuristic_score` | Implemented | Execution-tested | Not applicable | **D18 corrected 2026-11-02**: stored as `orf["heuristic_score"]`; `calculate_orf_confidence` remains as an alias of the former name. A ranking aid, not a probability. |
| — | ORF density, frame statistics | `calculate_orf_density`, `calculate_frame_statistics` | Implemented | Execution-tested | Not applicable | |
| — | ORF translation | `translate_orf_sequence`, `initiation_codons`, `is_initiation_codon` | Implemented | Execution-tested | **Independently validated** (Biopython, tables 1/2/4/11) | **D1 corrected 2026-10-31**: initiator Met only for initiation codons of the genetic code (table 1: ATG/CTG/GTG/TTG); internal codons use the elongation code. **D2 corrected 2026-10-31**: accepts the ORF dict and uses its `genetic_code`; a conflicting table raises; a bare string defaults to table 1. Agrees with `core.api` for every ORF found in tables 1/2/4/11. |
| — | Start-codon classification and scores | `classify_start_codon`, `calculate_cds_quality`, `calculate_orf_heuristic_score` | Implemented | Execution-tested | Not applicable | **S19 corrected 2026-10-31**: all three use the initiation codons of the ORF's genetic code, so a CTG ORF found with `table_starts` is classified and scored like an equivalent GTG ORF. |
| — | Circular / origin-spanning ORFs | — | **Not implemented** | Not applicable | Not applicable | `core.api.find_orfs` adds a limitation note for circular records. |

Test coverage: 19 test references — the best-covered module.

---

## 4. RNA analysis — `analyses/rna.py` (2,617 lines, 47 functions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| — | Source resolution (ORF vs whole sequence) | `resolve_rna_source` | Implemented | Integration-tested | Not applicable | **Confirmed intact**: `mode="orf"`/label "Selected ORF" with a selection, whole sequence and "Selected DNA" without. |
| — | Transcription | `transcribe_rna` | Implemented | Execution-tested | Not applicable | T→U on the given strand. |
| 1 | RNA Overview | `print_rna_overview`, `analyze_rna` | Implemented | Execution-tested | Not applicable | |
| 2 | Translation & Codon Analysis | `translate_rna_frame`, `analyze_rna_three_frames`, `analyze_rna_codon_usage` | Partially implemented | Execution-tested | Not applicable | **Always NCBI table 1**, even for a table-2/4/11 ORF (documented P1-9 limitation; `main.py` warns). Not part of D1/D2; the remaining place where a selected ORF is not translated with its own genetic code. |
| 3 | RNA Features | `find_rna_hairpins`, `summarize_rna_structure`, `find_rna_start_stop_landscape`, `find_rna_homopolymers` | Implemented | Execution-tested | Not applicable | Hairpins are a complementarity proxy, **not** a folding model. P2-5 G·U pairing and optional RNAfold deferred. |
| 4 | Complete RNA Report | `print_rna_report` | Implemented | Reachable | Not applicable | |
| 5 | RNA Sequence view | `print_rna_sequence_view` | Implemented | Reachable | Not applicable | |
| — | Base counts, GC/AU content, skew, profiles | `count_rna_bases`, `calculate_rna_gc_content`, `calculate_rna_gc_profile`, … | Implemented | Execution-tested | Not applicable | |
| — | Molecular weight, complexity | `calculate_rna_molecular_weight`, `calculate_rna_complexity` | Implemented | Execution-tested | Not applicable | |
| — | Ambiguity, quality | `analyze_rna_ambiguity`, `analyze_rna_quality` | Implemented | Execution-tested | Not applicable | |
| — | Kozak / Shine–Dalgarno on RNA | `analyze_rna_kozak_context`, `analyze_rna_shine_dalgarno` | Implemented | Execution-tested | Not applicable | Duplicates the ORF-module logic on RNA alphabet. |
| — | Feature score | `calculate_rna_feature_score` | Implemented | Execution-tested | Not applicable | Heuristic. |
| — | Pipeline storage | `_store_rna_in_pipeline`, `analyze_rna_from_pipeline` | Implemented | Integration-tested | Not applicable | |

Test coverage: 2 test references.

---

## 5. Protein analysis — `analyses/protein.py` (3,777 lines, 55 functions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| — | Translation | `translate_protein`, `analyze_translation`, `resolve_codon` | Implemented | Execution-tested | Independently validated | GTG-ORF initiator → M matches `Bio.Seq.translate(table=11, cds=True)`. **D1 corrected 2026-10-31**: with `initiator_start=True` only an initiation codon of the genetic code becomes Met (previously any resolved first codon); same rule as the ORF translator. |
| — | Unresolved-residue handling | `get_positional_protein`, `get_resolved_protein` | Implemented | Execution-tested | Not applicable | Stage 5: `?` keeps its position; windows containing it are skipped, not shifted. |
| 1 | Protein Overview | `_print_protein_overview` | Implemented | Reachable | Not applicable | |
| 2 | Translation & Sequence | `_print_translation_and_sequence` | Implemented | Reachable | Not applicable | |
| 3 | Amino Acid Composition | `analyze_protein`, `calculate_amino_acid_percentages` | Implemented | Execution-tested | Not applicable | |
| 4 | Physicochemical & Structural Properties | see below | Implemented | Execution-tested | **Independently validated** | |
| — | Molecular weight | `calculate_protein_molecular_weight` | Implemented | Execution-tested | Independently validated | Ubiquitin 8564.84 Da = ExPASy. Returns a dict with min/max/ambiguity flags. |
| — | Isoelectric point | `calculate_isoelectric_point` | Implemented | Execution-tested | Independently validated | 6.5616 vs ProtParam 6.56; lysozyme 9.3238 vs 9.32. |
| — | GRAVY | `calculate_gravy` | Implemented | Execution-tested | Independently validated | −0.48947 vs −0.489. |
| — | Instability index | `calculate_instability_index` | Implemented | Execution-tested | Independently validated | 36.055 vs 36.06; 16.091 vs 16.09. **IG weight corrected 2026-11-02** (−7.49 -> 1.0, the value in Biopython's table (`ProtParamData.DIWV`, which cites Guruprasad et al. 1990); the paper's table not consulted): the two audit proteins contain no IG dipeptide, so the error was invisible there; 20 of 31 random test proteins changed. Now within 2e-4 of Biopython on 31 proteins. |
| — | Aromaticity | `calculate_aromaticity` | Implemented | Execution-tested | Independently validated | Exact match. |
| — | Extinction coefficient | `calculate_extinction_coefficient` | Implemented | Execution-tested | Independently validated | Reduced/oxidized both match, incl. 4 disulfide pairs in lysozyme. |
| — | Charge at pH | `calculate_charge_at_ph`, `_get_pka_tables` | Implemented | Execution-tested | Not applicable | |
| — | Hydropathy profile, hydrophobic regions | `calculate_hydropathy_profile`, `find_hydrophobic_regions` | Implemented | Execution-tested | Not applicable | Kyte–Doolittle. |
| — | Transmembrane screen | `find_transmembrane_domains` | Implemented | Execution-tested | Not applicable | **A screen, not TMHMM/Phobius** — correctly labeled. |
| — | Signal-peptide screen | `predict_signal_peptide` | Implemented | Execution-tested | Not applicable | **Not SignalP** — correctly labeled. |
| — | Flexibility profile | `calculate_flexibility_profile` | Implemented | Execution-tested | Not applicable | |
| — | Secondary-structure propensity | `calculate_secondary_structure_propensity`, `…_fraction` | Implemented | Execution-tested | Not applicable | Propensity scales, not a predictor. |
| — | Aliphatic index | `calculate_aliphatic_index` | Implemented | Execution-tested | Not applicable | |
| 5 | Complete Protein Report | `_print_complete_protein_report` | Implemented | Reachable | Not applicable | |
| — | Quality assessment | `analyze_protein_quality`, `validate_protein_sequence` | Implemented | Execution-tested | Not applicable | |

Test coverage: 8 test references, incl. the 5 Stage-5 positional tests.

---

## 6. Restriction analysis — `analyses/restriction.py` (2,230 lines, 57 functions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| — | Enzyme database (22 enzymes) | `data/restriction_enzymes.py` | Implemented | Execution-tested | **Independently validated** | 22/22 match REBASE data via Biopython on site, both cut offsets, polarity, overhang length. |
| — | Site finding, CONFIRMED/POSSIBLE | `find_restriction_site_matches` | Implemented | Execution-tested | Independently validated | Supports `topology="circular"` incl. origin-spanning. |
| — | Legacy site list | `find_restriction_sites` | Implemented | Execution-tested | Not applicable | **D3 corrected 2026-10-31**: `topology` parameter. |
| 2 | Restriction Map | `create_restriction_map` | Implemented | Execution-tested | **Independently validated** | **D3 corrected 2026-10-31**: circular maps include origin-spanning sites (`wraps_origin`), coordinates modulo length; the menu asks for the topology. Circular map vs Biopython on pUC19: 132/132. |
| — | Cut positions | `calculate_cut_positions`, `get_cut_positions` | Implemented | Execution-tested | Independently validated | Validates `cut_top`/`cut_bottom` against the site length. |
| 3 | Digest & fragments (linear) | `perform_digest`, `calculate_linear_fragments` | Implemented | Execution-tested | **Independently validated** | pUC19, 22 enzymes: 0 differences vs Biopython. |
| 3 | Digest & fragments (circular) | `calculate_circular_fragments` | Implemented | Execution-tested | **Independently validated** | Including two origin-rotation stress cases: 0 differences. |
| 3 | Gel bands | `calculate_gel_bands` | Implemented | Execution-tested | Not applicable | Groups equal-size fragments. |
| 1 | Restriction Overview | `_restriction_overview_ui`, `analyze_restriction` | Implemented | Execution-tested | Not applicable | D3 corrected: map, per-enzyme results, unique/non-cutters, internal sites and pairs follow the selected topology and agree with the digest. |
| — | Unique cutters / non-cutters / internal sites | `find_unique_cutters`, `find_non_cutters`, `find_internal_sites` | Implemented | Execution-tested | Not applicable | |
| — | Overhangs and polarity | `get_overhang_sequence`, `_site_overhang_polarity` | Implemented | Execution-tested | Independently validated | Derived from cut coordinates, not a hand table. |
| — | Sticky-end compatibility | `sticky_ends_are_compatible`, `analyze_cloning_compatibility` | Implemented | Execution-tested | Independently validated | Polarity enforced (P0-3). |
| — | Enzyme-pair recommendation | `recommend_enzyme_pairs` | Implemented | Execution-tested | Not applicable | |
| — | Double digest | `analyze_double_digest`, `run_double_digest` | Implemented | Execution-tested | Not applicable | |
| — | Database validation | `validate_enzyme_database` | Implemented | Execution-tested | Not applicable | |
| — | Methylation notes | enzyme `methylation_note` | Implemented | Execution-tested | Not applicable | Display only; blocking and star activity not modeled. |
| 4 | Complete Restriction Report | `print_restriction_report` | Implemented | Reachable | Not applicable | |
| — | Topology selection UI | `_synapse_ui_select_topology` | Implemented | Reachable | Not applicable | |
| — | Type IIS enzymes | — | **Not implemented** | Not applicable | Not applicable | P2-7. No Type IIS enzyme is in the database, so no wrong result arises today; a reverse-strand search is required before one is added. |
| — | Self-test | `run_restriction_self_test` | Implemented | Execution-tested | Not applicable | |

---

## 7. Cloning — `analyses/cloning.py` (4,351 lines, 71 functions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| — | Site finding | `find_restriction_sites` → `restriction.find_recognition_site_matches` | Implemented | Execution-tested | Independently validated (shared matcher, pUC19 vs Biopython) | **D4 corrected 2026-10-31**: the restriction matcher is used; only CONFIRMED sites are returned (`include_possible=True` for POSSIBLE); the private `_recognition_matches` was removed. **D5 corrected 2026-10-31**: `topology` parameter. |
| — | Enzyme analysis, map | `analyze_enzyme`, `build_restriction_map` | Implemented | Execution-tested | Not applicable | D4 corrected: `possible_sites` / `possible_site_count` reported; an enzyme with POSSIBLE sites never counts as unique (ends, single-enzyme, suitability, planner diagnostics). D5 corrected: `topology`, recorded in the result. |
| — | Vector / insert topology | `vector_topology`, `insert_topology` on every cloning entry point | Implemented | **Integration-tested** | Not applicable | **D5 corrected 2026-10-31**: origin-spanning sites of circular vectors and inserts are found; digests use the molecule's topology; pipeline and menus pass the loaded record's topology. Assumption stated in results: on a circular molecule the segment between the two sites in coordinate order is replaced/excised. Tests: `CloningCircularTopologyTests`. |
| — | Overhang / polarity accessors | `get_overhang`, `get_overhang_polarity`, `get_cut_type`, `is_blunt_enzyme` | Implemented | Execution-tested | **Independently validated** | Derived from cut coordinates; 484-pair matrix matches Biopython. |
| — | End compatibility | `ends_compatible`, `analyze_end_compatibility`, `_end_signature`, `_site_polarity` | Implemented | Execution-tested | **Independently validated** | Returns `None` for undetermined polarity — a deliberate tri-state. |
| — | Junction analysis | `analyze_junction`, `_end_definition` | Implemented | Execution-tested | Not applicable | |
| 1 | Automatic Cloning Planner | `automatic_cloning_planner`, `find_best_cloning_strategies`, `_find_all_cloning_strategies` | Implemented | Integration-tested | Not applicable | D4 and D5 corrected 2026-10-31. **S18**: slow with all 22 enzymes (minutes on a 1.3 kb vector). |
| 2 | Single-Enzyme Analysis | `_cloning_ui_single_enzyme`, `_print_single_enzyme_report` | Implemented | Execution-tested | Not applicable | **Corrected 2026-11-03**: the report raised `KeyError: 'left'` for every enzyme (since at least the audit); found by `tests/test_rendered_menus.py`. |
| 3 | Two-Enzyme Strategy Analysis | `analyze_enzyme_pair`, `analyze_directional_strategy`, `_cloning_ui_two_enzyme_strategy` | Implemented | Execution-tested | Not applicable | Directional-cloning logic incl. internal-site checks. |
| — | Strategy scoring / ranking | `_score_strategy`, `score_cloning_strategy`, `rank_cloning_strategy`, `compare_cloning_strategies` | Implemented | Execution-tested | Not applicable | Heuristic. **D4 corrected 2026-10-31**: only CONFIRMED sites define cloning ends. |
| — | Digest fragments | `_digest_fragments`, `_selected_digest` | **Partially implemented** | Execution-tested | Not applicable | D5 corrected: sites and fragments follow the molecule's topology. **S17**: cut positions are taken as the recognition start, so fragment lengths are off by the cut offset (e.g. PstI by 4 nt). |
| — | Enzyme suitability / compatible enzymes | `assess_enzyme_suitability`, `find_compatible_enzymes` | Implemented | Execution-tested | Not applicable | |
| — | Diagnostics when no strategy exists | `_site_diagnostics`, `_no_strategy_diagnostics` | Implemented | Execution-tested | Not applicable | Good failure reporting. |
| — | Design reasons, ASCII orientation view | `get_design_reasons`, `ascii_orientation_visualization` | Implemented | Reachable | Not applicable | |
| — | Reports | `print_cloning_report`, `_print_single_strategy_report`, `print_planner_report` | Implemented | Reachable | Not applicable | |
| — | Self-tests | `run_cloning_self_test`, `run_cloning_ui_self_test` | Implemented | Execution-tested | Not applicable | |

---

## 8. Primer design — `analyses/primer.py` (5,412 lines, 126 definitions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| — | Nearest-neighbor Tm | `calculate_tm_nearest_neighbor` (one definition since **D16**), `calculate_duplex_thermodynamics` | Implemented | Execution-tested | **Independently validated** | **40/40 comparisons vs primer3-py in the audit, max difference 0.003 °C; extended in Stage I to 2,916 comparisons over Mg²⁺ 0-3 mM, max 0.0032 °C, kept in the suite.** |
| — | SantaLucia 1998 / 2004 / legacy parameter sets | `_nn_parameter_set`, `_nearest_neighbor_sum` | Implemented | Execution-tested | Independently validated | Default 1998 unified set. |
| — | Salt correction, Mg²⁺ model | `calculate_salt_correction`, `calculate_effective_ionic_strength`, `_effective_monovalent_salt_mM` | Implemented | Execution-tested | Independently validated | von Ahsen 2001 monovalent-equivalent; `na_only` and legacy 4× selectable. |
| — | Symmetry / self-complementarity term | `_is_self_complementary`, `calculate_duplex_concentration_term` | Implemented | Execution-tested | Independently validated | C/4 vs C and −1.4 e.u. handled. |
| — | Wallace Tm | `calculate_tm_wallace` | Implemented | Execution-tested | Not applicable | Offered as a UI alternative; correctly described as approximate. |
| — | Salt-adjusted Tm | `calculate_tm_salt_adjusted` | Implemented | Execution-tested | Not applicable | |
| — | Tm model comparison | `compare_tm_methods`, `get_thermodynamic_summary` | Implemented | Execution-tested | Not applicable | Good transparency feature. |
| 1 | Design Primers | `design_primers`, `design_primers_professional`, `generate_primer_candidates` | Implemented | Execution-tested | Not applicable | |
| — | Candidate scoring and ranking | `_v5`–`_v11` `*_pair_quality` / `*_rank_professional_primer_pairs` chain | Implemented | Execution-tested | Not applicable | **D16 corrected 2026-10-31**: the reachability review showed only the last binding of each name was reachable; superseded generations removed. Merging the surviving generations is not done. |
| — | GC clamp, homopolymer, repeat screens | `count_gc_clamp`, `longest_homopolymer`, `longest_tandem_repeat` | Implemented | Execution-tested | Not applicable | |
| — | Hairpin / self-dimer / heterodimer | `analyze_hairpin`, `analyze_self_dimer`, `analyze_heterodimer`, `_alignment_stem`, `_three_prime_complementarity` | Implemented | Execution-tested | **Not validated** | Heuristic complementarity screens, **not** thermodynamic structure prediction. Correctly disclaimed in METHODS. |
| — | Specificity screen | `screen_primer_specificity` | Implemented | Execution-tested | Not applicable | Against supplied sequence only, not a genome (P3-2). |
| — | Pairing and amplicon | `build_primer_pair`, `get_amplicon`, `summarize_primer_pair` | Implemented | Execution-tested | Not applicable | |
| 2 | Single Primer Analysis | `analyze_primer`, `analyze_primer_professional`, `_primer_ui_single_analysis` | Implemented | Execution-tested | Not applicable | |
| — | Failure diagnostics | `diagnose_thermodynamic_failure`, `_professional_failure_diagnostics` | Implemented | Execution-tested | Not applicable | |
| — | Condition entry UI | `_primer_ui_read_conditions`, `_primer_ui_select_tm_model` | Implemented | Reachable | Not applicable | |
| — | Restriction-site-driven primer design | — | **Not implemented** | Not applicable | Not applicable | Brief §10.2 describes the workflow; no wiring from restriction-site selection to primer design exists. |
| — | Reports | `print_primer_report`, `print_professional_primer_report`, `get_professional_primer_summary` | Implemented | Reachable | Not applicable | |
| — | Self-tests | `run_primer_self_test`, `run_thermodynamic_self_test`, `run_primer_v3_self_test`, `run_primer_v3_self_test_final`, `run_primer_v4_screening_self_test` | Implemented | Execution-tested | Not applicable | Several self-test entry points, a remnant of the generations removed under **D16**. |

---

## 9. Mutation analysis — `analyses/mutation.py` (4,184 lines, 82 definitions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| 1 | Compare Sequences | `compare_sequences`, `global_align`, `build_alignment_blocks` | Implemented | Execution-tested | Not applicable | Needleman–Wunsch, match 2 / mismatch −1 / gap −2. **S9 quantified 2026-10-31**: O(n·m) memory, ~45 bytes per cell (about 4.5 GB for 10 kb x 10 kb); the note states it. |
| — | Mutation detection | `detect_mutations`, `detect_mutations_from_alignment`, `_detect_simple_indel` | Implemented | Execution-tested | Not applicable | Pure length changes use a prefix/suffix representation. |
| — | Substitution classification | `classify_substitution`, `is_transition` | Implemented | Execution-tested | Not applicable | |
| — | Coding consequences | `analyze_coding_consequences`, `classify_coding_substitution`, `annotate_mutations_with_coding_consequences` | Implemented | Execution-tested | Not applicable | Synonymous / missense / nonsense / start-loss / stop-loss / frameshift / in-frame indel all covered by self-tests. |
| — | Start-loss restricted to codon 1 | `classify_coding_substitution(codon_index=)` | Implemented | Execution-tested | Not applicable | P0-6 holds. Limitation: alternative initiation codons in the reference are not recognized as starts. |
| — | CDS from the selected ORF | `orf_reference_cds`, `assess_coding_frame` | Implemented | Integration-tested | Not applicable | Both strands; variant CDS derived only at equal length. |
| 2 | Mutation Statistics | `calculate_mutation_statistics` | Implemented | Execution-tested | Not applicable | |
| — | Mutation spectrum, hotspots | `build_mutation_spectrum`, `find_mutation_hotspots` | Implemented | Execution-tested | Not applicable | |
| — | HGVS-*like* notation | `mutation_hgvs_c`, `mutation_hgvs_p`, `_mutation_notation` | Partially implemented | Execution-tested | **Not validated** | Explicitly labeled HGVS-like; not checked against den Dunnen et al. 2016. |
| — | Indel left-normalization | `normalize_indel_left` | Implemented | Execution-tested | Not applicable | Separate from the variant module's implementation — duplication. |
| — | Context extraction | `extract_mutation_context` | Implemented | Execution-tested | Not applicable | |
| — | Restriction-site impact | `restriction_site_impact` | Implemented | Execution-tested | Not applicable | Cross-module check. |
| — | Primer / CRISPR overlap | `primer_mutation_overlap`, `crispr_mutation_overlap` | Implemented | Execution-tested | Not applicable | P2-4's cross-check, partially realized. |
| — | Protein comparison | `compare_proteins` | Implemented | Execution-tested | Not applicable | |
| 3 | Coding Consequence Analysis (UI) | `run_mutation_analysis_ui` branch 3 | Implemented | Reachable | Not applicable | |
| 4 | Mutation Report | `build_mutation_report`, `print_mutation_report`, `advanced_mutation_report`, `mutation_map` | Implemented | Reachable | Not applicable | |
| — | Advanced analysis | `create_advanced_analysis`, `mutation_signature` | Implemented | Execution-tested | Not applicable | |
| — | JSON export | `result_to_json` | Implemented | Execution-tested | Not applicable | |
| — | Validation | `validate_mutation_result` | Implemented | Execution-tested | Not applicable | |
| 5 | Self-test | `run_mutation_self_test` | Implemented | Execution-tested | Not applicable | 44/44 incl. the TEST_8 family. |
| — | Size guard | 2,000 nt note | Partially implemented | Execution-tested | Not applicable | **S9**: a note with the measured figures, not a limit. P2-6 10 kb target unmet. |

---

## 10. Variant analysis — `analyses/variant_analysis.py` (2,322 lines, 73 definitions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| 1 | Load VCF | `parse_vcf`, `parse_vcf_file`, `_load_vcf_text_interactive` | Implemented | Execution-tested | Not applicable | VCF 4.x text, gzip and BGZF. |
| — | Header `Number=` parsing | `parse_vcf_header_numbers`, `_field_number` | Implemented | Execution-tested | Not applicable | VCF 4.3 reserved keys as fallback. |
| — | Multi-allelic splitting | `_select_allele_specific_info`, `_remap_sample_data_for_alt`, `_genotype_index` | Implemented | Execution-tested | Not applicable | A/R/G cardinality; G uses F(j,k) = k(k+1)/2 + j. Other ALT alleles in GT → 0, as `bcftools norm -m-`. |
| — | Genotype parsing, phasing, categories | `parse_genotype`, `genotype_category` | Implemented | Execution-tested | Not applicable | |
| — | FILTER semantics | `_normalize_filter_status` | Implemented | Execution-tested | Not applicable | `.` = not applied, not failed. |
| 2 | Add Variant Manually | `_prompt_variant` | Implemented | Reachable | Not applicable | |
| 3 | Normalize Variants | `left_normalize_variant`, `normalize_variants`, `trim_common_bases` | Implemented | Execution-tested | **Partially validated** | Algorithm hand-verified against Tan et al. 2015 on 5 constructed cases (audit §3.7). **No bcftools comparison** — blocked on Windows. |
| — | REF validation | `check_reference_alleles`, `_reference_allele_matches` | Implemented | Execution-tested | Not applicable | **D12 corrected 2026-10-31**: a contig without a reference is flagged `REFERENCE_CONTIG_NOT_FOUND` and listed in `contigs_not_found`; `not_checked` is split by reason. Previously it was skipped silently. Tests: `tests/test_phase6_tier1.py::VariantContigMatchingTests`. |
| — | Contig-name matching | `resolve_contig_names`, `ContigResolution`, `contig_matching=` | Implemented | Execution-tested | Not applicable | Exact by default; opt-in `chr_prefix` pairs `1`/`chr1` (flag `CONTIG_NAME_ALIASED`, CHROM unchanged). Mitochondrial names never paired; case variants not folded. **Not covered**: contig-length checks against `##contig`, RefSeq/GenBank accession names. |
| 4 | Filter Variants | `variant_passes_filter`, `filter_variants`, `VariantFilter` | Implemented | Execution-tested | Not applicable | Missing-value filtering covered by a self-test. |
| 5 | Variant Statistics | `calculate_variant_statistics`, `is_transition_variant` | Implemented | Execution-tested | Not applicable | Ti/Tv, allele frequencies, quality. |
| 6 | Analyze Variants | `analyze_variants`, `annotate_variant`, `annotate_from_info` | Implemented | Execution-tested | Not applicable | |
| 7 | Compare Variant Sets | `compare_variant_sets`, `variant_key_set` | Implemented | Execution-tested | Not applicable | |
| 8 | Cohort / Sample Analysis | `sample_variant_summary`, `calculate_cohort_summary` | Implemented | Execution-tested | Not applicable | |
| 9 | Export VCF | `variants_to_vcf`, `write_vcf` | Implemented | Execution-tested | Not applicable | Round-trip self-tested. |
| 10 | Export TSV | `variants_to_tsv` | Implemented | Execution-tested | Not applicable | |
| 11 | Full Variant Report | `build_variant_report`, `print_variant_report` | Implemented | Reachable | Not applicable | |
| — | Prioritization, intervals, grouping | `prioritize_variants`, `variants_in_interval`, `interval_variant_density`, `variants_by_*` | Implemented | Execution-tested | Not applicable | |
| — | Evidence container | `VariantEvidence`, `attach_evidence` | Implemented | Execution-tested | Not applicable | |
| — | Mutation-module bridge | `annotate_from_mutation_analysis` | Implemented | Execution-tested | Not applicable | |
| 12 | Self-test | `run_variant_analysis_self_test` | Implemented | Execution-tested | Not applicable | 27/27. |

---

## 11. CRISPR analysis engine — `analyses/crispr.py` (2,574 lines, 134 definitions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| 1 | Sequence & Region | `build_target_region`, `_load_orf_or_region`, `TargetRegion` | Implemented | Execution-tested | Not applicable | Region offsets verified: a sub-region at `region_start=61` reports parent-consistent coordinates. |
| — | ORF integration | `orf_target_context`, `extract_orf_context`, `core.orf_context.resolve_orf_region` | Implemented | **Integration-tested** | Not applicable | Minus-strand ORF → forward region 5–100 identical to the parent slice, strand `-`; 0 coordinate mismatches. P0-2 holds. |
| 2 | Batch Analysis | `batch_design_guides`, `batch_summary`, `load_fasta`, `parse_fasta_text` | Implemented | Execution-tested | Not applicable | |
| 3 | CRISPR System selection | `get_crispr_system`, `list_crispr_systems`, `format_system_details` | Implemented | Execution-tested | Not applicable | 4 systems: SpCas9 NGG/20, SpCas9-NG NG/20, SaCas9 NNGRRT/21, Cas12a TTTV/23. |
| 4 | Target Discovery | `find_targets`, `_find_targets_in_region`, `_candidate_from_oriented_window`, `pam_matches` | Implemented | **Execution-tested** | Not applicable | Both strands; 0 coordinate mismatches over all candidates in 3 sequences. Cas12a 5′-PAM geometry verified (PAM 41–44, guide 45–67). |
| — | Guide characterization | `characterize_guide`, `GuideMetrics`, `quality_flags` | Implemented | Execution-tested | Not applicable | GC, homopolymer, repeats, complexity, self-complementarity, hairpin proxy, ambiguity. **S4 resolved 2026-10-31**: thresholds sourced or labeled as Synapse conventions. |
| 5 | Guide Analysis & Ranking | `rank_guides`, `_component_*`, `ranking_mode_label` | Partially implemented | Execution-tested | Not applicable | Transparent weighted heuristic, correctly disclaimed. **D9 corrected 2026-10-31**: the on-target site no longer counts against specificity. |
| — | Cut-site / cleavage position | `cut_site`, `cut_site_complementary`, `cut_type` per guide | Implemented | Execution-tested | Not applicable | **D8 corrected 2026-10-31**: modeled cut per system with its citation, in the model, display and exports. |
| — | PAM quality | `pam_quality` | Implemented | Execution-tested | Not applicable | Pattern degeneracy only, explicitly not an activity score. |
| 6 | Mutation / Variant Targeting | `annotate_mutation_targets`, `mutation_position_in_guide`, `parse_mutation_target`, `parse_vcf_variant_line` | Implemented | Execution-tested | Not applicable | |
| 7 | Base Editing | `analyze_base_editing`, `analyze_base_editing_candidate`, `base_editing_window`, `get_base_editor` | Partially implemented | Execution-tested | Not applicable | 4 editor models (CBE/ABE/CBE+/ABE+). Guide-relative window verified (position 5 in window 4–8). **D10 corrected 2026-10-31**: a configured window is used. Bystander detection works. |
| 8 | Prime Editing | `design_prime_editing_candidates`, `_prime_edit_for_candidate` | Partially implemented | Execution-tested | Not applicable | **S3 corrected 2026-10-31**: the RTT is checked to contain the intended edit. Honest heuristic labeling. |
| 9 | Off-Target Analysis | `search_local_off_targets`, `search_local_off_targets_detailed`, `annotate_local_off_targets`, `_oriented_candidate_windows` | Partially implemented | **Execution-tested** | Not applicable | PAM-aware, seed-aware, mismatch positions reported. **D9 corrected 2026-10-31**: the on-target site is excluded. Scope limited to supplied sequences — correctly stated. |
| 10 | Specificity Profile | `specificity_profile`, `has_specificity_analysis` | Partially implemented | Execution-tested | Not applicable | **D9 corrected 2026-10-31**: `exact_matches` no longer counts the guide's own site. |
| — | **CFD scoring** | — | **Not implemented** | Not applicable | Not applicable | P2-3: deliberately rejected pending a cited published mismatch matrix. **Correct decision — do not fabricate weights.** |
| — | U6 poly-T rule | `quality_flags` poly-T warning | Partially implemented | Execution-tested | Not applicable | **S5 addressed 2026-10-31**: the flag is named as the Pol III (U6) terminator consideration and sourced; still a sequence rule, not a promoter model. |
| 11 | Reports & Export | `format_report`, `print_report`, `export_json`, `export_csv`, `export_guides_fasta`, `export_batch_summary`, `resolve_export_path` | Implemented | Execution-tested | Not applicable | **D11 corrected 2026-10-31**: exports go to `exports/crispr` in the project unless `SYNAPSE_EXPORT_DIR` is set. |
| — | Serialization | `guide_to_dict`, `result_to_dict`, `result_to_json` | Implemented | Execution-tested | Not applicable | Includes the modeled cut position (**D8 corrected 2026-10-31**). |
| 12 | Self-Test | `run_self_test` | Implemented | Execution-tested | Not applicable | 28/28. |
| 13 | Supported Systems | `list_crispr_systems` | Implemented | Execution-tested | Not applicable | |
| — | Config fields `include_low_quality`, `rank_by_gc`, `editing_window_start/end`, `base_editor` | `DesignConfig` | Implemented | Execution-tested | Not applicable | **D10 corrected 2026-10-31**: each field now acts, with the previous behavior as default (`CrisprConfigurationTests`). |
| — | Input validation | `normalize_sequence`, `validate_sequence` | Implemented | Execution-tested | Not applicable | Empty, invalid-symbol and all-N inputs all raise clear errors; short sequences warn and return 0 candidates. |

Test coverage: **4 test references for 2,574 lines** — the weakest coverage
relative to size; the 28 self-tests are the module's own fixtures.

---

## 12. BLAST search — `analyses/blast.py` (4,650 lines, 103 definitions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| 1 | Run BLAST Search (local) | `run_blast_analysis(engine="local")`, `LocalBlastEngine.run` | **Partially implemented** | **Execution-tested** | **Independently validated** | **Real BLAST+ 2.17.0+ run executed** (audit §6): 2 hits, coverage 100 % / 98.04 %, correct E-values and bands. The programmatic `local_database_path` branch now verifies the database prefix like a managed database (**D6**, **S7** corrected 2026-10-31; tested with real BLAST+). Protein queries (blastp/tblastn) were rejected by a type-vocabulary mismatch and the menu offered them no program (**D14-F2**, corrected 2026-10-31; a real local blastp run is now a test). |
| 1 | Run BLAST Search (remote) | `NCBIBlastClient.submit/status/fetch_xml/run` | Implemented | **Not tested** | Blocked by external dependency | No live call made in this audit. Privacy confirmation before sending is implemented. |
| — | Query building | `build_query`, `build_orf_query`, `infer_sequence_type`, `normalize_sequence` | Implemented | Execution-tested | Not applicable | ORF query retains global coordinates and strand. The query type comes from `core.alphabets`; an undetermined type is asked for (menu) or refused (API), RNA is searched as DNA with a note (**D14-F1**, corrected 2026-10-31). |
| — | FASTA query loading | `load_fasta_queries`, `_resolve_synapse_user_path` | Implemented | Execution-tested | Not applicable | |
| 2 | Resume Search by Request ID (RID) | `NCBIBlastClient.resume`, `core.api.blast_resume` | Implemented | Execution-tested (stubbed client) | Not applicable | **S6 implemented 2026-10-31** (P1-7). A live resume was not run (external request). |
| 3 | Search Configuration | `BlastParameters`, `.normalized()` | Implemented | Reachable | Not applicable | Program/database compatibility self-tested; program/query-type compatibility via `check_program_query_type` (D14-F2). |
| 4 | BLAST Database Manager | `blast_database_manager_menu`, `LocalBlastDatabaseManager` | Implemented | Execution-tested | Not applicable | 8 sub-entries: status, download/update, verify, create custom from FASTA, remove, information, BLAST+ tool status, change root. |
| 4.1/4.3/4.6 | Database status / scan / verify | `scan_database`, `verify_database`, `database_info`, `_parse_blastdbcmd_info` | Implemented | Execution-tested | Not applicable | 6 catalog entries correctly report `NOT_INSTALLED`. Multi-volume and alias discovery self-tested. |
| 4.2 | Download / update NCBI database | `update_ncbi_database` | Implemented | **Not tested** | Blocked by external dependency | Would be a very large download; deliberately not run. |
| 4.4 | Create custom database from FASTA | `create_custom_database` | Implemented | Not tested | Not applicable | `makeblastdb` is available, so this is now testable. |
| — | XML parsing | `parse_blast_xml`, `_parse_hsp`, `_parse_hit`, `_parse_parameters`, `_parse_database_statistics` | Implemented | **Execution-tested** | **Independently validated** | Validated on real BLAST+ output, not only synthetic XML. |
| 5 | Result Summary | `print_hit_summary`, `rank_hits` | Implemented | Execution-tested | Not applicable | |
| 6 | Hit Detail & Alignment | `print_hit_detail`, `alignment_text` | Implemented | Execution-tested | Not applicable | |
| — | Query coverage | `BlastHSP.query_coverage`, `BlastHit.best_coverage_percent` | Implemented | Execution-tested | **Independently validated** | Union of HSP query intervals / query length (NCBI "Query cover"); 51/51 and 50/51 verified. P0-1 confirmed on real output. |
| 7 | Evidence Profile | `calculate_evidence_score`, `build_evidence_profile`, `HitEvidence` | Implemented | Execution-tested | Not applicable | Bands gated at E ≤ 1e-3 (`DEFAULT_EVIDENCE_EVALUE_THRESHOLD`, configurable); `NOT_SIGNIFICANT` otherwise. |
| 8 | Taxonomic Overview | `group_hits_by_taxonomy`, `_taxonomy_source`, `_infer_scientific_name_from_definition` | Implemented | Execution-tested | Not applicable | Definition-line fallback is self-tested. |
| 9 | Compare BLAST Runs | `compare_blast_runs`, `BlastComparison` | Implemented | Execution-tested | Not applicable | |
| 10 | Export Results | `export_json`, `export_csv`, `export_text`, `run_to_dict`, `run_to_json` | Implemented | Execution-tested | Not applicable | |
| 11 | Full BLAST Report | `build_text_report` | Implemented | **Execution-tested** | Not applicable | Rendered from the real run. |
| — | Provenance | `BlastRunMetadata`, `database_details`, `reported_parameters` | Partially implemented | Execution-tested | Not applicable | Effective parameters, BLAST version, engine, local path and timestamps all recorded. Title/date/size are also recorded on the user-supplied-path branch (**S7** corrected 2026-10-31). |
| — | BLAST+ executable management | `BlastExecutableManager` (`available`, `resolve`, `require`, …) | Implemented | **Execution-tested** | Not applicable | Correctly located BLAST+ 2.17.0+ in `C:\Program Files\NCBI\`. |
| — | Short-query rescue | `_should_run_short_query_rescue`, task configuration | Implemented | Execution-tested | Not applicable | Self-tested policy for short queries. |
| 12 | Self-Test | `run_self_test` | Implemented | Execution-tested | Not applicable | 33/33. |
| 13/14 | Supported programs; database status | `PROGRAMS`, `LOCAL_DATABASE_CATALOG` | Implemented | Execution-tested | Not applicable | |
| — | Stored run artifact | `blast_results/TEST_8_REFERENCE_blast.json` | Legacy record | Not applicable | Not applicable | Labeled as a pre-P0-1 legacy record in `blast_results/README.md` (**D7** corrected 2026-10-31; deletion left to the owner). Originally: stale pre-P0-1 values (100 % coverage at E = 7.45, null query length) and an obsolete key schema; no test reads it. |

---

## 13. Literature search — `analyses/literature.py` (1,624 lines, 44 definitions)

| Menu | Feature | Function | Implementation | Execution | Validation | Notes |
|---|---|---|---|---|---|---|
| 1 | New Literature Search | `run_literature_search`, `LiteratureQuery` | Implemented | Execution-tested (offline) | Blocked by external dependency | No live call in this audit. |
| — | PubMed provider | `search_pubmed`, `_journal_from_pubmed_xml`, `_publication_types_from_pubmed_xml` | Implemented | Execution-tested (fixtures) | Not validated | Throttled 3/s, 10/s with `SYNAPSE_NCBI_API_KEY`. |
| — | Europe PMC provider | `search_europe_pmc` | Implemented | Execution-tested (fixtures) | Not validated | |
| — | OpenAlex provider | `search_openalex` | Implemented | Execution-tested (fixtures) | Not validated | Stage 3 fixed a malformed `YYYY-MM-01-01` date. |
| — | Crossref provider | `search_crossref` | Implemented | Execution-tested (fixtures) | Not validated | |
| — | Throttling and retry | `_throttle`, `_retry_delay`, `_request_text`, `_request_json` | Implemented | Execution-tested | Not applicable | Up to 3 attempts on 429/5xx, honours `Retry-After`; a simulated 503 is tested. |
| — | Date handling | `_date_bound`, `_date_in_range`, `_apply_date_validation`, `_future_date_note` | Implemented | **Execution-tested** | Not applicable | `2020-03` → 2020-03-01 / 2020-03-31; `2020` → 2020-12-31. Correct. |
| — | Deduplication | `deduplicate_articles`, `_is_placeholder_journal` | Implemented | Execution-tested (fixtures) | Not applicable | **S11**: cross-source duplicates unverified live (needs the owner's approval for external requests). |
| — | Author-name handling | `_normalize_person_name`, `_deduplicate_author_names` | Implemented | Execution-tested | Not applicable | |
| — | Relevance ranking | `rank_articles`, `_coverage`, `_tokenize`, `_query_phrase` | Implemented | Execution-tested | Not applicable | Text matching only — correctly stated not to measure evidence quality. |
| — | Metadata enrichment | `_enrich_missing_metadata` | Implemented | Execution-tested (fixtures) | Not applicable | |
| 2 | Search Configuration | `_interactive_menu` | Implemented | Reachable | Not applicable | Parameters entered per search. |
| 3 | Run Self-Test | `self_test` | Implemented | Execution-tested | Not applicable | Passes offline. |
| 4 | Data Sources & Scientific Limitations | `_interactive_menu` branch 4 | Implemented | Reachable | Not applicable | Good transparency feature. |
| — | Export | `export_results`, `article_to_dict` | Implemented | Execution-tested | Not applicable | |
| — | Sequence-context suggestion | `build_sequence_literature_context` | Implemented | Integration-tested | Not applicable | Offers the header as a suggestion; **no sequence data is sent**. |

---

## 14. Core — `core/` (8 modules, 1,806 lines)

| Component | Function/class | Implementation | Execution | Notes |
|---|---|---|---|---|
| Shared data model | `models.Interval`, `Feature`, `SequenceRecord`, `AnalysisResult` | Implemented | Execution-tested | One internal convention (0-based half-open + strand) with explicit 1-based conversion. `SequenceRecord.alphabet` comes from `core.alphabets.classify_sequence` with an evidence level (**D14 corrected 2026-10-31**; level renamed from "confidence" under D18). |
| Hashing / provenance | `models.sha256_text`, `sha256_file`, `provenance()` | Implemented | Execution-tested | Present in every stored result. |
| Non-interactive API | `api.load_records`, `sequence_overview`, `find_orfs`, `analyze_rna`, `analyze_protein_from_cds`, `analyze_protein`, `restriction_digest`, `design_primers`, `compare_sequences`, `analyze_vcf`, `design_crispr_guides` | **Partially implemented** | Execution-tested | 16 entry points, no `input()`. Cloning, variant filtering, mutation reporting, BLAST and literature added 2026-10-31 (**S14**); the two that make external requests refuse to run without an explicit confirmation argument. Input types checked with `core.alphabets`, optional `sequence_type=` (**D14-F7** corrected 2026-10-31). **Missing**: cloning, variant filtering, mutation reporting, BLAST, literature (**S14**). `analyze_vcf` now reports `normalized: True` only when every variant had a reference, with `data["normalization"]` and per-contig warnings (**D12** corrected 2026-10-31). |
| Pipeline context and store | `pipeline.PipelineContext.store`, `store_result`, `latest`, `results` | Implemented | **Integration-tested** | ORF change clears `rna`, `protein`, `crispr`, `blast` — verified. `completed` has a fixed shape and `run_orf` takes the ORF search settings (**D19** corrected 2026-10-31). |
| ORF-dependent invalidation | `pipeline.ORF_DEPENDENT_MODULES` | Implemented | Integration-tested | Confirmed working. |
| RNA / protein source resolution | `pipeline.get_rna_source_dna`, `get_protein_source_rna` | Implemented | **Integration-tested** | Selected-ORF-vs-whole-sequence rule confirmed intact. |
| Pipeline module runners | `pipeline.run_sequence`, `run_orf`, `run_rna`, `run_protein`, `run_restriction`, `run_cloning`, `run_cloning_planner`, `run_core_pipeline` | Implemented | Integration-tested | **D19 corrected 2026-10-31**: `run_orf` takes the menu's search settings. Workflow tests (10.1, 10.7) run these runners end to end. |
| ORF hand-over adapter | `orf_context.resolve_orf_region`, `OrfRegion`, `is_synapse_orf` | Implemented | **Integration-tested** | Validates the ORF against the parent and refuses parent coordinates on mismatch. P0-2's central fix. |
| Sequence alphabets and type classifier | `alphabets.classify_sequence`, `require_sequence_type`, `normalize_sequence`, `declared_type_from_molecule` | Implemented | Execution-tested | Letter sets **independently validated** against `Bio.Data.IUPACData`. Added in Phase 6 (D14/D13). Used by `fileio`, `SequenceRecord` (D14-F3, F11) and BLAST (D14-F1). |
| Shared sequence tools | `seqtools.normalize_dna`, `complement`, `reverse_complement`, `gc_fraction`, `GeneticCode`, `codon_dict`, `get_genetic_code` | Implemented | Execution-tested | **Independently validated** against `Bio.Data.CodonTable`. `normalize_dna` has the same contract as `core.validation` (**D13** corrected 2026-10-31). |
| Legacy validation | `validation.normalize_dna`, `validate_dna`, `get_dna_errors` | Implemented | Execution-tested | Thin layer over `core.alphabets`; empty input invalid (**D13**, **D14-F4** corrected 2026-10-31). |
| Reproducible report | `report.build_analysis_report`, `report_to_markdown`, `export_analysis_report`, `to_jsonable`, `canonical_json` | Implemented | **Execution-tested** | JSON and Markdown both written; per-analysis payload SHA-256; rejects other suffixes. **S8 corrected 2026-10-31**. |
| Logging | `logging_setup.configure_logging` | Implemented | Execution-tested | Rotating `logs/synapse.log`; `main.report_error` writes tracebacks there. |

---

## 15. Support modules

| Component | Function | Implementation | Execution | Notes |
|---|---|---|---|---|
| FASTA loading | `fileio.fasta` | Implemented | Execution-tested | `topology=circular` header token supported. |
| GenBank loading | `fileio.genbank` | **Partially implemented** | Execution-tested | LOCUS topology read (circular fixture confirmed, 3 features). Locations: simple spans, `complement()`, `join()`; **S10 corrected 2026-10-31**: other forms keep the raw string and the record says so in a note. |
| Loader dispatch | `fileio.loader.load_sequences` | Implemented | Execution-tested | Content-sniffs; clear errors for empty/header-less files; records path, format, SHA-256, index. **D14 corrected 2026-10-31**: declared molecule type read and every record classified. |
| Main menu rendering | `ui.menu.show_menu`, `select_gene` | Implemented | Execution-tested | **D17 corrected 2026-11-01**: 70-column rules from `ui.style`; `[0] Exit`. |
| Sequence overview rendering | `ui.report.print_results` (728 lines) | Implemented | Execution-tested | **D17 corrected 2026-11-01**: one width through `ui.style`. |
| Enzyme data | `data.restriction_enzymes` (22 enzymes + 6 accessors) | Implemented | Execution-tested | **Independently validated** (22/22). `ENZYME_DATABASE_VERSION = "synapse-2026-10-29"`. |
| Amino-acid data | `data.amino_acids` | Implemented | Execution-tested | |
| Genetic code data | `data.genetic_code` | Implemented | Execution-tested | Derived from `core.seqtools.codon_dict` since Stage 4. |

---

## 16. Advertised-but-absent features

As found by the audit, with the outcome of Phase 6 in each row.

| Advertised where | Feature | Status |
|---|---|---|
| `main.py` menu `[11]` | AI Experiment Planner | **Not implemented**, and since 2026-10-31 (**D15**) the menu says so instead of "coming soon". P3-1. |
| `crispr.DesignConfig` | `include_low_quality`, `rank_by_gc`, `editing_window_start`, `editing_window_end`, `base_editor` | Implemented 2026-10-31 (**D10**): each field acts (`CrisprConfigurationTests`). |
| CRISPR menu `[4]`/`[5]` | Cut position for the chosen system | Implemented 2026-10-31 (**D8**), with the source per system. |
| CRISPR menu `[9]`/`[10]` | On-target exclusion from off-target counts | Implemented 2026-10-31 (**D9**). |
| `restriction` docs / P2-7 | Type IIS enzymes, reverse-strand site search | **Not implemented** — no Type IIS enzyme in the database, so no wrong output today. |
| `variant` P1-5 | Contig-name unification (`chr1` vs `1`) | **Implemented 2026-10-31** (**D12**): opt-in `contig_matching="chr_prefix"`; unmatched contigs are flagged, not skipped. |
| BLAST P1-7 | RID resume | Implemented 2026-10-31: `resume()`, menu [2] (numbered [14] until the Stage F renumbering), `core.api.blast_resume()`. Tested against a stub; a live resume was not run. |
| ORF METHODS | 5′-partial ORFs; ORFs across a circular origin | **Not implemented** (documented). |
| Brief §10.2 | Restriction-site-driven primer design | **Not implemented** as a feature: no menu passes a chosen site to primer design. The two modules share one coordinate frame (integration test 10.2). |
| P2-3 | CFD off-target scoring | **Not implemented** — deliberately deferred pending a citable matrix. Correct call. |

**No original feature was lost.** Every capability listed in
`docs/AUDIT_STATUS.md` and in the previous cycle's reports is accounted for
above, either as working, as defective with a `D`/`S` reference, or as
explicitly not implemented.

---

## 17. Validation ledger

### Performed in this audit (external references)

| # | Validation | Reference | Result |
|---|---|---|---|
| 1 | 22 restriction enzymes | Biopython 1.88 `Bio.Restriction` (REBASE/EMBOSS) | 22/22 on 5 properties |
| 2 | Digests on pUC19 (L09137.2) | Biopython, 22 enzymes × 2 topologies × 3 rotations | 132 comparisons, 0 differences |
| 3 | End-compatibility matrix | Biopython `ovhg`/`ovhgseq` | 484 pairs, 0 mismatches |
| 4 | Melting temperature | primer3-py 2.4.0 (SantaLucia) | 40 comparisons, max 0.003 °C |
| 5 | Genetic codes 1/2/4/11 | `Bio.Data.CodonTable` | Identical, incl. start sets |
| 6 | Protein properties | `Bio.SeqUtils.ProtParam`, ExPASy | pI, GRAVY, instability, aromaticity, extinction, MW all match |
| 7 | Local BLAST pipeline | BLAST+ 2.17.0+ on the bundled database | Real run; coverage, E-values, bands, provenance correct |
| 8 | Variant left-alignment | Tan et al. 2015, hand-verified | 5/5 cases correct, incl. correct non-shifts |
| 9 | Circular restriction map (after D3) | Biopython 1.88 on pUC19, 22 enzymes × 2 topologies × 3 rotations | 132 comparisons, 0 differences |
| 10 | Melting temperature, extended (Stage I, 2026-11-02) | primer3-py 2.4.0, 81 primers × 36 conditions (Mg²⁺ 0-3 mM) | 2,916 comparisons, max 0.0032 °C; in the suite |
| 11 | Enzymes, digests, end compatibility (Stage I) | Biopython 1.88 `Bio.Restriction`, 6-kb test sequence, 2 topologies × 3 rotations | 22/22 enzymes, 132 digests and 484 pairs, 0 differences; in the suite |
| 12 | Protein descriptors (Stage I) | Biopython 1.88 `ProtParam`, 31 proteins | found the IG instability weight error (corrected); then all within tolerance; in the suite |

### Still open

| Validation | Blocker | Route |
|---|---|---|
| `bcftools norm` comparison | No native Windows build | WSL or a Linux host |
| ORFfinder comparison | No local binary | NCBI ORFfinder standalone, or a curated RefSeq CDS set |
| NEBcutter comparison | Web service | Substantively superseded by validation 2 |
| Live remote NCBI BLAST | Not run by choice | One confirmed run recording RID and provenance |
| Local BLAST on a real database | `nt` is ~100 GB+ | `swissprot` (~0.5 GB) as a proportionate step |
| Python 3.10 suite | 3.10 not installed (only 3.14.6); installing an interpreter needs the owner's approval | Install 3.10 and run the suite. Meanwhile `pyright --pythonversion 3.10` (112 diagnostics, baseline 134) and `Python310CompatibilityTests` (3.10 grammar, pre-3.12 f-string rules) run; the latter found and fixed one CRISPR import failure on 3.10/3.11 |
| HGVS conformance | No validator | `hgvs` package or Mutalyzer; currently labeled "HGVS-like" |
| Primer hairpin/dimer screens | Heuristic by design | `primer3.calc_hairpin` / `calc_heterodimer` are available and would be a direct comparison |

---

## 18. How to keep this document honest

1. Never move a row to **Independently validated** without naming the
   external reference and the observed numbers.
2. **Execution-tested** requires that the output was inspected, not merely
   that no exception was raised.
3. A green test suite does not change any status in this table. All 19
   confirmed defects live in code the 109 tests report as passing.
4. When a defect is fixed, cite the test that would now catch it.
5. Re-check the "advertised-but-absent" table (§16) at the end of Phase 6;
   every row must be implemented, removed from the interface, or
   explicitly labeled unavailable to the user.
