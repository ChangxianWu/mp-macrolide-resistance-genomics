# Data directory

Public *Mycoplasma pneumoniae* assemblies are **not** redistributed in this
repository. Download them from NCBI (or ENA) using the accession list in
Supplementary Table S1 of the manuscript.

The first pipeline step (`scripts/01_build_manifest.py`) produces a download
manifest that maps ENA accessions to assembly accessions. Place the downloaded
FASTA files in this directory using the naming convention expected by
`scripts/02_download_fasta.py`, or simply run that script to fetch them.

Everything in this directory except this README is git-ignored.
