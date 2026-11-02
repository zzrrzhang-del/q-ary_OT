# Reproducing the q-ary OT Finite-Batch Artifacts

Run all commands from the repository root. See the top-level `README.md` for
the tested dependency versions and the recommended run order.

## Validate representative model results

```powershell
python .\scripts\validate_examples.py
```

Expected values for `delta=8`, `kappa=s=128`:

- `q=2`: `n_r=255`, `n_r*r=255`
- `q=4`: `n_r=170`, `n_r*r=340`
- `q=16`: `n_r=136`, `n_r*r=544`
- `q=256`: `n_r=128`, `n_r*r=1024`

The validation also checks finite-batch break-even behavior:
with `B_base=100` bits, `q=4` does not beat the binary baseline for `delta=8`.

## Run the code search

The paper uses cached selected code candidates from:

- `tex\generated\data\code_search_all.csv`
- `tex\generated\data\code_search_selected.csv`
- `tex\generated\data\code_search_metadata.json`
- `tex\generated\code_search_table.tex`

To search local constructions and query CodeTables where the public table
covers the parameter range, run:

```powershell
python .\scripts\code_search.py --force --online-codetables
```

For an offline search using local construction families only, run:

```powershell
python .\scripts\code_search.py --force
```

The offline search combines juxtaposed simplex, shortened BCH parameters from
q-cyclotomic coset degree accounting, repetition codes for `k=1`, and
projective-line Griesmer codes for `k=2`. For non-binary BCH, the scripts record
cyclotomic-coset degree and designed-distance guarantees rather than emitting
generator polynomials or matrices in Python. Magma verifies selected codes.

## Run Magma construction certificates

Magma can verify the selected `delta=12` candidates. A fresh clone retrieves
the GF(4) generator matrix from CodeTables during this step; see the top-level
`THIRD_PARTY_DATA.md` for its provenance and retrieval requirements.

```powershell
python .\scripts\run_magma_certificates.py --magma "D:\Program\Magma\Magma\magma.exe"
```

This constructs and checks:

- the shortened binary BCH candidate `[313,12,135]_2`;
- the CodeTables `[177,6,128]_4` generator matrix;
- the shortened non-binary BCH candidates `[184,4,135]_8` and `[153,3,135]_16`;
- the explicit projective-line `[135,2,132]_64` construction;
- the repetition `[132,1,132]_{4096}` construction.

Generated Magma certificate artifacts:

- `tex\generated\magma_certificates_table.tex`
- `tex\generated\data\magma_certificates.csv`
- `tex\generated\data\magma_certificates.log`
- `tex\generated\data\magma_certificates_input.m` (generated locally, not tracked)
- `tex\generated\data\codetables_q4_177_6_matrix.csv` (downloaded locally, not tracked)

The Magma certificate step verifies the constructed selected candidates; it
does not perform a global best-known-code search. The `q=4` candidate is
obtained from CodeTables and checked locally using its stored generator matrix.

## Model implemented by the scripts

- Candidate radices are all `r | delta`.
- The statistical parameter is padded as `s_r = r * ceil(s / r)` so that the
  CCG subfield-check condition `r | s_r` holds.
- The distance target is `d_star = max(kappa, s_r)`.
- The selected code length is the shortest available candidate among the local
  construction/search families and cached CodeTables witnesses.
- The Griesmer length bound is converted into a communication lower bound:

```text
B_lower,r(M) = g_r * (B_base + r * (M + 4s_r)) + 2s_r M + 2s_r delta
```

- Workload points are classified as certified binary, certified non-binary
  advantage, or code-dependent by comparing selected witnesses with these lower
  bounds.
- The finite-batch communication model is:

```text
B_ext,r(M) = (M + 4s_r)n_r r + 2s_r M + 2s_r delta
B_total,r(M) = n_r B_base + B_ext,r(M)
```

## Generate paper tables and figures

```powershell
python .\scripts\generate_artifacts.py
```

Generated manuscript tables:

- `tex\generated\construction_witness_table.tex`
- `tex\generated\magma_certificates_table.tex`
- `tex\generated\critical_thresholds_table.tex`
- `tex\generated\padding_sensitivity_table.tex`
- `tex\generated\application_cases_table.tex`
- `tex\generated\data\params_all.csv`
- `tex\generated\data\construction_witnesses.csv`
- `tex\generated\data\code_search_all.csv`
- `tex\generated\data\code_search_selected.csv`
- `tex\generated\data\magma_certificates.csv`
- `tex\generated\data\critical_thresholds.csv`
- `tex\generated\data\certification_points.csv`
- `tex\generated\data\base_ot_points.csv`
- `tex\generated\data\breakeven_all.csv`
- `tex\generated\data\breakeven_summary.csv`
- `tex\generated\data\application_cases.csv`
- `tex\generated\data\base_ot_sensitivity.csv`
- `tex\generated\data\padding_sensitivity.csv`

Generated figures:

- `tex\figures\pareto_delta4.pdf`
- `tex\figures\pareto_delta8.pdf`
- `tex\figures\pareto_delta12.pdf`
- `tex\figures\pareto_overview.pdf`
- `tex\figures\radix_regions_delta4.pdf`
- `tex\figures\radix_regions_delta8.pdf`
- `tex\figures\radix_regions_delta12.pdf`
- `tex\figures\radix_regions_overview.pdf`
- `tex\figures\analysis_framework.pdf`
- `tex\figures\base_ot_sensitivity.pdf`

The `radix_regions_delta*.pdf` files are certified-region phase diagrams, not
plain selected-radix heatmaps. Their sample CSV files retain both the
certificate class and the selected witness radix.

The scripts use deterministic integer formulas and fixed logarithmic grids. No
random seed is needed. `B_base` values are sensitivity points anchored by
reported single-base-OT communication where available.

## Manuscript build

This repository contains the generated tables and figures, but not the full
manuscript source or its LLNCS class file. Compile the paper in the separate
manuscript project after copying the generated artifacts there.
