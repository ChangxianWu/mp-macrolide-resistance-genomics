# Publishing this repository on GitHub (maintainer notes)

These steps are for the corresponding author. Nothing in this folder has been
uploaded anywhere yet. A descriptive repository name such as
`mp-macrolide-resistance-genomics` is recommended.

## 1. Pre-flight check

From the repository root:

- **No data**: `data/` contains only `README.md`; `.gitignore` excludes everything
  else in it. No assemblies, result tables, figures, DOCX/PDF/ZIP files.
- **No secrets**: no tokens, no local absolute paths (scripts resolve the root
  relative to `scripts/`).
- **No patient information**: all genomes are public NCBI assemblies.
- **Expected content**: `README.md`, `LICENSE`, `requirements.txt`,
  `CITATION.cff`, `.gitignore`, `run_pipeline.py`, `docs/`, `data/README.md`, and
  35 `.py` files under `scripts/`.

Optional check:

```bash
git init -b main
git add .
git status          # review the file list before committing
```

## 2. Create and push the repository

Create an empty repository on GitHub (no README, no license), then:

```bash
git init -b main
git add .
git commit -m "Genomic epidemiology of macrolide resistance in public Mycoplasma pneumoniae assemblies"
git remote add origin https://github.com/<account>/<repository>.git
git push -u origin main
```

If GitHub asks for credentials, use a personal access token or an SSH key — do
not paste a token into any project file.

## 3. After the push

On GitHub, check that:

- the README renders (title, layout tree, stage table, code blocks);
- the file tree contains `LICENSE`, `requirements.txt`, `CITATION.cff`,
  `run_pipeline.py`, `scripts/` (35 scripts), `docs/`, `data/README.md`;
- no data, results, figures or documents were uploaded;
- optional: fill in the About description and topics
  (`mycoplasma-pneumoniae`, `macrolide-resistance`, `genomics`, `23s-rrna`), and
  create a `v1.0.0` release.

## 4. Link the repository in the manuscript

After the push, record the repository URL in one place in the manuscript project
(a single constant in the manuscript configuration, outside this repository) and
regenerate the manuscript variants. The Data availability statement then points
to this repository while the per-step output tables remain in the supplementary
data set. Keep the repository and the frozen submission archive consistent.

## 5. License

The repository is released under the MIT License (`LICENSE`, copyright 2026 the
five authors). The same license is declared in `CITATION.cff`; keep the two files
in sync if the license ever changes.

## 6. Consistency

The supplementary Data Set S1 zip is the frozen, submission-time archive (code +
per-step output tables + MANIFEST). The GitHub repository is the living code
copy. If the code changes after submission, update the repository and note the
change in the manuscript's revision letter.
