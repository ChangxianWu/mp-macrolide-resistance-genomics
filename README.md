# Genomic epidemiology of macrolide resistance in public *Mycoplasma pneumoniae* assemblies, 1944–2025

Analysis code accompanying the manuscript:

> **A 25-SNP Genomic-Similarity Group Carries Most Macrolide Resistance–Associated
> Genotypes in Public *Mycoplasma pneumoniae* Genomes, 1944–2025**
> Changxian Wu, Min Wang, Fan Liu, Xiaoman Zheng, Hong Su
> (manuscript in preparation)

This is a dry-lab re-analysis of 646 public genome assemblies. The repository
contains the scripted pipeline that produces every number, table and main-text
figure in the manuscript. No sequencing data, intermediate caches or result
tables are committed.

---

## Repository layout

```
.
├── README.md
├── LICENSE                 MIT license
├── requirements.txt        pinned Python dependencies
├── CITATION.cff            machine-readable citation metadata
├── run_pipeline.py         single entry point (stage order, dry-run, subset runs)
├── scripts/                36 pipeline scripts (run order = numeric prefix)
├── docs/PIPELINE.md        stage-by-stage commands, outputs and method notes
└── data/README.md          where to place downloaded assemblies (git-ignored)
```

The flat `scripts/` folder is deliberate: every script resolves the project root
as the parent of `scripts/`, so numbered scripts can be run in order without
installation. `run_pipeline.py` groups them into nine stages.

**Not included**: raw assemblies, caches, generated tables, figures and
documents; superseded analysis scripts; the DOCX typesetting scripts. Result
tables can be obtained from the corresponding author on reasonable request.

---

## Data

All genomes are **public NCBI assemblies**; the complete per-genome accession
list, metadata and typing results are in **Supplementary Table S1** of the
manuscript. Sparse records were completed from linked ENA BioSample entries
(114 Japanese assemblies, 1976–2024).

Run `scripts/01_build_manifest.py` first to generate the download manifest, then
either place the assemblies in `data/` or let `scripts/02_download_fasta.py`
download them. Everything under `data/` except `README.md` is git-ignored.

---

## Install and run

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate    |    macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt

python run_pipeline.py --list          # nine stages, 36 scripts
python run_pipeline.py --dry-run       # print the exact commands
python run_pipeline.py --from 03_reference   # resume from a stage
python run_pipeline.py --stage 09_figures_tables
```

Each script can also be run directly, e.g. `python scripts/14_coresnp.py`.

---

## Pipeline stages

| # | Stage | Scripts | What it produces |
|---|---|---|---|
| 1 | `01_panel_metadata` | `01`, `01b` | NCBI/ENA download manifests |
| 2 | `02_acquisition` | `02` (`02b` alternative) | assemblies in `data/` |
| 3 | `03_reference` | `03e` | calibrated 23S anchor coordinates |
| 4 | `04_typing` | `04`, `06`, `07`, `08`, `09`, `10` | 23S, P1 and MLST typing, assembly QC, integrated matrix, frequency tables |
| 5 | `05_coresnp_tree` | `11`, `12`, `14`, `15` | P1 clusters, k-mer NJ tree, anchor-chain core-SNP matrix (742,308 positions), SNP tree |
| 6 | `06_temporal` | `15b`, `16`, `17`, `18`, `19`, `19b`, `20`, `21` | clock diagnostics, TreeTime exports and five-setting grid, root-free Fitch reconstruction, descriptive temporal models |
| 7 | `07_group_structure` | `23`, `26`, `32`, `33`, `35` | 25-SNP single-linkage group (352 assemblies), full-tree ancestry, NNI/pruning sensitivity, within-group substructure, genome-jackknife stability |
| 8 | `08_normalization` | `28` | display-layer region names (`Mainland China`, `Taiwan, China`); statistics unchanged |
| 9 | `09_figures_tables` | `22`, `24`, `25`, `29`, `30`, `31`, `34` | main figures, manuscript tables, Supplementary Table S1 |

Run stage 8 before stage 9 so that labels and tables use the normalized display
layer (the archived analysis was produced in that order).

---

## Method notes

- **23S typing**: dual 22-nt anchors flanking positions 2063/2064/2617, both
  strands, calls only when all anchor pairs agree. These are
  resistance-associated genotypes; no MIC or culture susceptibility was measured.
- **Anchor-chain core-SNP calling**: unique 25-mer index of the M129 reference,
  copy-number-aware matching in both orientations, same-strand chaining,
  Needleman–Wunsch for unequal gaps, indel neighbourhoods masked.
- **Two invariants**: anchor bodies always carry the *reference* base; every
  k-mer classifier is canonical (`min(forward, reverse complement)`).
- **Interpretation limits**: single-linkage groups are similarity units, not
  clades or transmission chains; Fitch change directions follow the displayed
  arbitrary root; TreeTime rates were not identifiable, so no molecular clock or
  absolute dates are reported.

---

## Citation, contact, license

Please cite the accompanying manuscript (`CITATION.cff` is provided; the journal
reference will be added after publication).

Corresponding author: Changxian Wu, College of Agriculture and Animal Husbandry
Science and Technology, Ningxia Vocational and Technical College for Minorities,
Wuzhong 751100, China — wchangx173@163.com.

The code is released under the MIT License (see `LICENSE`).
