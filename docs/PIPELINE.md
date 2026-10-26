# Pipeline run order

Run everything from the repository root. `run_pipeline.py` is the entry point;
each script can also be executed directly. Outputs are written under `results/`
and `manuscript/`, which are git-ignored.

```bash
python run_pipeline.py --list
python run_pipeline.py --dry-run
python run_pipeline.py                 # all nine stages
python run_pipeline.py --from 04_typing
```

Every script resolves the repository root as the parent of `scripts/`.

---

## 1. Panel metadata

```bash
python scripts/01_build_manifest.py     # NCBI assembly_name -> download manifest
python scripts/01b_ena_manifest.py      # ENA Portal cross-check
```

## 2. Data acquisition

```bash
python scripts/02_download_fasta.py --included-only --workers 6
# alternative mirror:
# python scripts/02b_download_ena.py --included-only --workers 5
```

Downloads are idempotent and resumable. Assemblies are not redistributed.

## 3. Reference calibration

```bash
python scripts/03e_calibrate_sites.py
```

Produces `reference/reference_anchors.json` and the calibration report. An
earlier EC-numbering route was superseded and is not included.

## 4. Typing, QC and the integrated matrix

```bash
python scripts/04_type_23S.py        # positions 2063/2064/2617
python scripts/06_summary_tables.py  # frequency by year/country
python scripts/07_w2qc.py            # N50/GC and sketch redundancy
python scripts/08_p1_type.py         # P1 type (M129 / FH)
python scripts/09_mlst.py            # eight-locus MLST
python scripts/10_matrix.py          # integrated per-genome matrix
```

Main outputs: the per-genome typing tables and the integrated matrix under
`results/`, plus frequency tables by year and country.

## 5. Core-genome SNP matrix and tree

```bash
python scripts/11_p1_cluster.py
python scripts/12_mash_nj.py
python scripts/14_coresnp.py         # anchor-chain calling (~40 min)
python scripts/15_snp_matrix.py      # 742,308 core positions; NJ tree
```

Main outputs: per-genome allele vectors, pairwise SNP distances and the NJ tree
under `results/coresnp/`.

## 6. Temporal and ancestral analyses

```bash
python scripts/15b_clock.py
python scripts/16_export_aln.py
python scripts/17_datetree.py
python scripts/18_export_full.py
python scripts/19_treetime.py
python scripts/19b_treetime_report.py
python scripts/20_resistance_origin.py
python scripts/21_epi_dynamics.py
```

TreeTime is run under five rooting/coalescent-prior settings. The rate estimates
are not identifiable; no substitution rate or absolute divergence time is
reported. Fitch reconstruction is root-free; the directions of individual changes
follow the displayed root.

## 7. Genomic-similarity group and sensitivity analyses

```bash
python scripts/23_epidemic_clade.py
python scripts/26_fulltree_ancestry.py
python scripts/32_tree_perturbation_fitch.py
python scripts/33_cluster_substructure.py
```

Produces the 25-SNP single-linkage group (352 assemblies), within-group
substructure, NNI topology perturbation and panel-pruning sensitivity.

## 8. Display-layer normalization

```bash
python scripts/28_country_norm.py
```

Rewrites region labels to `Mainland China` / `Taiwan, China` / `United States` /
`Not specified` in the frequency table and the manuscript table CSV. Statistics
do not change; the NCBI raw country string remains in Supplementary Table S1.

## 9. Figures, tables and supplementary table

```bash
python scripts/22_figures_en.py
python scripts/24_manuscript_tables.py
python scripts/25_supplementary.py
python scripts/29_fig3_redesign.py
python scripts/30_figS1_threshold.py
python scripts/31_figS2_composition.py
python scripts/34_fig5_temporal_series.py
```

Outputs: main and supplementary figures under `manuscript/figures/`, manuscript
tables under `manuscript/tables/`, and the per-genome Supplementary Table S1.

---

## Invariants that must not be changed

1. **Anchor bodies carry the reference base.** Filling reverse anchors with the
   complement produced ~50% spurious differences in reverse-oriented contigs.
2. **All k-mer classifiers are canonical** (`min(forward, reverse complement)`).
   Non-canonical counts missed reverse-oriented contigs and produced false
   "divergent" P1 lineages.
3. **Single-linkage groups are similarity units.** They are not monophyletic
   groups, transmission chains, or evidence of a single ancestral mutation.
4. **No molecular clock.** Pairwise distance–year slopes are descriptive; the
   TreeTime grid was not identifiable.

## Excluded from this repository

Superseded analysis scripts — the first EC-numbering calibration route, the
first-version working figures, the review-response patch and the bibliography
utility — together with the DOCX typesetting scripts are retained in the frozen
Supplementary Data Set S1 archive.
