# Finite-Batch q-ary OT Extension: Analysis Artifacts

This repository contains the analysis code and fixed data used for the paper
*Finite-Batch q-ary OT Extension for High-Arity Private Lookup: Certified
Radix--Code Co-Design*. It reproduces the finite-batch communication model,
code-candidate selection, certificate-region calculations, tables, and figures.
It is **not** an implementation or benchmark of the CCG OT protocol or a
complete private-lookup application.

## Contents

| Path | Purpose |
| --- | --- |
| `scripts/radix_cost_model.py` | Security padding, code parameters, Griesmer bounds, communication costs, and break-even thresholds. |
| `scripts/code_search.py` | Simplex, BCH, repetition, and projective-line candidates; optional CodeTables queries. |
| `scripts/run_magma_certificates.py` | Magma construction and minimum-distance checks for selected `delta=12` codes. |
| `scripts/generate_artifacts.py` | Generates the CSV data, LaTeX tables, and PDF figures. |
| `scripts/validate_examples.py` | Checks the paper's representative parameter and finite-batch examples. |
| `scripts/REPRODUCIBILITY.md` | Detailed formulas, artifact list, and optional verification steps. |
| `tex/generated/data/` | Frozen scalar search inputs, provenance metadata, Magma transcript, and result CSVs. |
| `tex/generated/*.tex`, `tex/figures/*.pdf` | Generated manuscript tables and figures. |
| `THIRD_PARTY_DATA.md` | CodeTables attribution, matrix retrieval, and redistribution scope. |

The manuscript source and third-party reference PDFs are not part of this
artifact package. The generated tables and figures can be copied into a
separate LaTeX paper project.

## Versions and dependencies

- Tested on Python **3.13.5**, NumPy **2.1.3**, and Matplotlib **3.10.0**.
  Install the pinned Python dependencies from `requirements.txt`.
- The main model, search, validation, and figure generation use Python and its
  standard library plus NumPy and Matplotlib. No network access is needed for
  the main reproduction path below.
- **Optional:** Magma **V2.12-9** is used to rerun the selected code-distance
  checks. Install Magma separately and pass the executable path with `--magma`.
  The runner uses PowerShell on Windows and standard input on other systems.
  A fresh clone also needs network access to fetch the GF(4) matrix from
  CodeTables; that matrix is not included here.
- **Optional:** `code_search.py --online-codetables` and
  `run_magma_certificates.py` query CodeTables and require network access.
  The frozen scalar data records CodeTables access on 2026-08-10.

## Reproduce the paper artifacts

Run commands from the repository root. Create a Python 3.13 virtual
environment with `python -m venv .venv`, activate it for your shell, and then
run:

```bash
python -m pip install -r requirements.txt
python scripts/validate_examples.py
python scripts/generate_artifacts.py
```

For example, on Windows PowerShell, activate the environment with
`.\.venv\Scripts\Activate.ps1`; on macOS/Linux use
`source .venv/bin/activate`. The generator reads the frozen
`tex/generated/data/code_search_selected.csv` and
`tex/generated/data/code_search_metadata.json`, then writes tables and CSVs
under `tex/generated/` and figures under `tex/figures/`. Keeping the frozen
candidate files is necessary to reproduce the paper's selected-code results
without querying a live best-known-code table.

## Optional verification and new searches

To rerun the six selected code-distance checks on Windows, install Magma and
pass its executable path:

```bash
python scripts/run_magma_certificates.py --magma "D:\Program\Magma\Magma\magma.exe"
```

On another system, replace the path with the local Magma executable. The script
downloads the GF(4) matrix from CodeTables on a fresh clone and writes a local
matrix CSV, `magma_certificates_input.m`, `magma_certificates.log`, and
`magma_certificates.csv`. The matrix CSV and Magma input are Git-ignored. The
included transcript and CSV let readers inspect the recorded checks without
installing Magma. See `THIRD_PARTY_DATA.md` for source attribution, the pinned
matrix digest, and the optional TLS workaround.

Fresh searches are optional and can change the selected-code CSVs. Run them in
a separate checkout if you want to retain the frozen paper results:

```bash
python scripts/code_search.py --force
python scripts/code_search.py --force --online-codetables
```

The first command searches local construction families only. The second also
queries the current public CodeTables pages, so its results can differ from the
frozen 2026-08-10 inputs. A fresh search should be followed by
`python scripts/generate_artifacts.py`.

## Scope and provenance

The reported comparisons concern communication within the CCG q-ary OT
extension family. CodeTables entries provide selected code candidates;
Griesmer bounds provide the lower-bound certificates. The Magma checks verify
the selected `delta=12` constructions, not global code-length optimality.
The base-OT communication values are sensitivity anchors, not measurements of
an end-to-end protocol implementation.

This repository's original code, documentation, and generated artifacts are
released under the [MIT license](LICENSE), copyright (c) 2026 Zhejiang Lab.
The third-party CodeTables matrix is not included or relicensed; its source,
original construction credit, access date, and retrieval dependency are
documented in `THIRD_PARTY_DATA.md`. Cite CodeTables alongside this repository.
When citing the artifacts in the paper, use an immutable release or commit
reference.
