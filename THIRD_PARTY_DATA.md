# CodeTables data provenance

The selected GF(4) code of length 177 and dimension 6 uses a generator matrix
published on [CodeTables](https://www.codetables.de/BKLC/BKLC.php?k=6&n=177&q=4).
That page attributes the construction to **Axel Kohnert** and is maintained by
**Markus Grassl**. The fixed search metadata records access on **2026-08-10**.
CodeTables provides [citation guidance](https://www.codetables.de/citation.php).

We found no explicit matrix redistribution license on the construction or
citation pages when checked on 2026-11-03. The matrix is therefore **not
distributed in this repository** and is **not covered by the Zhejiang Lab MIT
license**. The generated Magma input file also embeds the matrix and is not
distributed. The included `magma_certificates.csv` and
`magma_certificates.log` record verification results, not matrix entries.

To rerun the GF(4) distance check, `scripts/run_magma_certificates.py` retrieves
the matrix from the CodeTables URL above, validates its 6-by-177 shape and field
symbols, and checks its canonical SHA-256 digest against
`12252741684441389e88b10aa8c0dfd34858ae9916a6b105aebba349d91727cf`.
The matrix CSV and Magma input are generated locally and ignored by Git. This
step requires network access to CodeTables and a working Magma installation.
If the page becomes unavailable or its matrix changes, the GF(4) verification
cannot be reproduced from this repository alone. The main cost-model, table,
and figure reproduction uses the frozen scalar candidate data and stays offline.

The fetch uses certificate-verified HTTPS by default. If the local TLS trust
store rejects the CodeTables certificate, the optional
`--allow-unverified-tls` flag skips certificate verification for this request;
the pinned SHA-256 check still rejects a different matrix. Use this option
only when certificate-verified retrieval is unavailable.
