"""Search q-ary code candidates for the OT parameter model.

The search combines:

1. Juxtaposed simplex codes.
2. Shortened BCH candidates, using q-cyclotomic cosets to compute the generator
   degree and BCH designed-distance guarantee.
3. Exact Griesmer constructions for k=1 and k=2:
   repetition codes for k=1 and projective-line Griesmer codes for k=2.
4. Optional CodeTables best-known-table probes where the public table
   covers the target field size and length range.

For non-binary prime-power fields, the BCH search records cyclotomic-coset
degree and designed-distance guarantees rather than constructing generator
polynomials in Python. Magma verifies the selected codes independently. Online
CodeTables probes are cached as scalar CSV artifacts and are not required for
the offline BCH/simplex search.

Use `run_magma_certificates.py` to construct and verify the selected delta=12
code candidates with the local Magma executable.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable

from radix_cost_model import RadixParams, divisors, griesmer_bound, params_for


ROOT = Path(__file__).resolve().parents[1]
TEX_DIR = ROOT / "tex"
GENERATED_DIR = TEX_DIR / "generated"
DATA_DIR = GENERATED_DIR / "data"

CODETABLES_URL = "https://www.codetables.de/BKLC/BKLC.php?k={k}&n={n}&q={q}"
CODETABLES_MAX_N = {
    2: 256,
    3: 243,
    4: 256,
    5: 130,
    7: 100,
    8: 130,
    9: 130,
    16: 130,
}
BOUNDS_RE = re.compile(
    r"lower bound:</TD><TD>\s*([0-9]+)\s*</TD>.*?"
    r"upper bound:</TD><TD>\s*([0-9]+)\s*</TD>",
    re.IGNORECASE | re.DOTALL,
)
METHOD_PRIORITY = {
    "repetition": 0,
    "projective-line Griesmer": 0,
    "juxtaposed simplex": 1,
    "shortened BCH": 2,
    "CodeTables witness": 3,
}


@dataclass(frozen=True)
class CodeCandidate:
    delta: int
    kappa: int
    s_target: int
    r: int
    q: int
    k: int
    s_r: int
    d_star: int
    method: str
    n_r: int
    d_actual: int
    bit_proxy: int
    griesmer_bound: int
    griesmer_gap: float
    distance_overshoot: float
    parent_n: int | None = None
    parent_k: int | None = None
    bch_start: int | None = None
    bch_degree: int | None = None
    external_lower_bound: int | None = None
    external_upper_bound: int | None = None
    source_url: str = ""
    verified: str = ""
    note: str = ""


def tex_escape(value: object) -> str:
    return (
        str(value)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("&", r"\&")
        .replace("%", r"\%")
    )


def candidate_from_repeated_simplex(p: RadixParams) -> CodeCandidate:
    return CodeCandidate(
        delta=p.delta,
        kappa=p.kappa,
        s_target=p.s_target,
        r=p.r,
        q=p.q,
        k=p.k,
        s_r=p.s_r,
        d_star=p.d_star,
        method="juxtaposed simplex",
        n_r=p.n_r,
        d_actual=p.d_actual,
        bit_proxy=p.bit_proxy,
        griesmer_bound=p.griesmer_bound,
        griesmer_gap=p.griesmer_gap,
        distance_overshoot=p.distance_overshoot,
        parent_n=p.L_r,
        parent_k=p.k,
        verified="formula",
        note="explicit construction",
    )


def candidate_from_repetition(p: RadixParams) -> CodeCandidate | None:
    if p.k != 1:
        return None
    n = p.d_star
    g = griesmer_bound(p.q, p.k, p.d_star)
    return CodeCandidate(
        delta=p.delta,
        kappa=p.kappa,
        s_target=p.s_target,
        r=p.r,
        q=p.q,
        k=p.k,
        s_r=p.s_r,
        d_star=p.d_star,
        method="repetition",
        n_r=n,
        d_actual=p.d_star,
        bit_proxy=n * p.r,
        griesmer_bound=g,
        griesmer_gap=n / g,
        distance_overshoot=1.0,
        verified="explicit",
        note="optimal for k=1",
    )


def candidate_from_projective_line(p: RadixParams) -> CodeCandidate | None:
    if p.k != 2:
        return None
    n = griesmer_bound(p.q, p.k, p.d_star)
    return CodeCandidate(
        delta=p.delta,
        kappa=p.kappa,
        s_target=p.s_target,
        r=p.r,
        q=p.q,
        k=p.k,
        s_r=p.s_r,
        d_star=p.d_star,
        method="projective-line Griesmer",
        n_r=n,
        d_actual=p.d_star,
        bit_proxy=n * p.r,
        griesmer_bound=n,
        griesmer_gap=1.0,
        distance_overshoot=1.0,
        verified="theorem",
        note="Griesmer optimal for k=2",
    )


def cyclotomic_coset_ids(n: int, q: int) -> tuple[list[list[int]], list[int]]:
    seen = [False] * n
    cosets: list[list[int]] = []
    ids = [-1] * n
    for a in range(n):
        if seen[a]:
            continue
        x = a
        coset: list[int] = []
        while not seen[x]:
            seen[x] = True
            ids[x] = len(cosets)
            coset.append(x)
            x = (x * q) % n
        cosets.append(coset)
    return cosets, ids


def window_degrees(ids: list[int], coset_sizes: list[int], width: int) -> Iterable[tuple[int, int]]:
    n = len(ids)
    counts = [0] * len(coset_sizes)
    degree = 0
    for offset in range(width):
        coset_id = ids[offset % n]
        if counts[coset_id] == 0:
            degree += coset_sizes[coset_id]
        counts[coset_id] += 1
    for start in range(n):
        yield start, degree
        outgoing = ids[start % n]
        counts[outgoing] -= 1
        if counts[outgoing] == 0:
            degree -= coset_sizes[outgoing]
        incoming = ids[(start + width) % n]
        if counts[incoming] == 0:
            degree += coset_sizes[incoming]
        counts[incoming] += 1


def bch_coset_best(
    p: RadixParams,
    max_parent_n: int | None = None,
    designed_distance: int | None = None,
) -> CodeCandidate | None:
    """Find the shortest shortened BCH candidate by coset-degree accounting."""

    d = designed_distance or p.d_star
    max_n = max_parent_n or p.n_r
    best: CodeCandidate | None = None
    for parent_n in range(d, max_n + 1):
        if math.gcd(parent_n, p.q) != 1:
            continue
        cosets, ids = cyclotomic_coset_ids(parent_n, p.q)
        coset_sizes = [len(c) for c in cosets]
        for start, degree in window_degrees(ids, coset_sizes, d - 1):
            parent_k = parent_n - degree
            if parent_k < p.k:
                continue
            shortened_n = parent_n - (parent_k - p.k)
            g = griesmer_bound(p.q, p.k, p.d_star)
            candidate = CodeCandidate(
                delta=p.delta,
                kappa=p.kappa,
                s_target=p.s_target,
                r=p.r,
                q=p.q,
                k=p.k,
                s_r=p.s_r,
                d_star=p.d_star,
                method="shortened BCH",
                n_r=shortened_n,
                d_actual=d,
                bit_proxy=shortened_n * p.r,
                griesmer_bound=g,
                griesmer_gap=shortened_n / g,
                distance_overshoot=d / p.d_star,
                parent_n=parent_n,
                parent_k=parent_k,
                bch_start=start,
                bch_degree=degree,
                verified="coset-degree",
                note="BCH designed-distance guarantee",
            )
            if is_better_candidate(candidate, best):
                best = candidate
                if candidate.n_r == candidate.griesmer_bound and candidate.distance_overshoot == 1:
                    return best
    return best


def query_codetables(q: int, k: int, n: int, timeout: float) -> tuple[int, int, str] | None:
    url = CODETABLES_URL.format(q=q, k=k, n=n)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            html = response.read().decode("utf-8", errors="replace")
    except (OSError, TimeoutError, urllib.error.URLError):
        return None
    match = BOUNDS_RE.search(html)
    if match is None:
        return None
    lower, upper = (int(match.group(1)), int(match.group(2)))
    return lower, upper, url


def codetables_best(
    p: RadixParams,
    timeout: float = 5.0,
    end_n: int | None = None,
) -> CodeCandidate | None:
    """Find the shortest CodeTables witness in the covered public table range."""

    max_table_n = CODETABLES_MAX_N.get(p.q)
    if max_table_n is None:
        return None
    start = p.griesmer_bound
    end = min(p.n_r, max_table_n, end_n or p.n_r)
    if start > end:
        return None
    for n in range(start, end + 1):
        bounds = query_codetables(p.q, p.k, n, timeout)
        if bounds is None:
            continue
        lower, upper, url = bounds
        if lower >= p.d_star:
            return CodeCandidate(
                delta=p.delta,
                kappa=p.kappa,
                s_target=p.s_target,
                r=p.r,
                q=p.q,
                k=p.k,
                s_r=p.s_r,
                d_star=p.d_star,
                method="CodeTables witness",
                n_r=n,
                d_actual=lower,
                bit_proxy=n * p.r,
                griesmer_bound=p.griesmer_bound,
                griesmer_gap=n / p.griesmer_bound,
                distance_overshoot=lower / p.d_star,
                external_lower_bound=lower,
                external_upper_bound=upper,
                source_url=url,
                verified="online table",
                note=f"lower={lower}; upper={upper}; accessed={date.today().isoformat()}",
            )
    return None


def candidate_key(candidate: CodeCandidate) -> tuple[int, float, int, str]:
    return (
        candidate.n_r,
        candidate.distance_overshoot,
        METHOD_PRIORITY.get(candidate.method, 9),
        candidate.method,
    )


def is_better_candidate(candidate: CodeCandidate, incumbent: CodeCandidate | None) -> bool:
    if incumbent is None:
        return True
    return candidate_key(candidate) < candidate_key(incumbent)


def all_candidates_for(
    delta: int,
    kappa: int,
    s_target: int,
    r: int,
    use_codetables: bool = False,
    codetables_timeout: float = 15.0,
) -> list[CodeCandidate]:
    p = params_for(delta, kappa, s_target, r)
    candidates = [candidate_from_repeated_simplex(p)]
    for optional in [candidate_from_repetition(p), candidate_from_projective_line(p)]:
        if optional is not None:
            candidates.append(optional)
    bch = bch_coset_best(p)
    if bch is not None:
        candidates.append(bch)
    if use_codetables:
        local_best = min(candidates, key=candidate_key)
        codetable = None
        if local_best.n_r > p.griesmer_bound:
            codetable = codetables_best(
                p,
                timeout=codetables_timeout,
                end_n=local_best.n_r - 1,
            )
        if codetable is not None:
            candidates.append(codetable)
    return candidates


def selected_candidates(
    deltas: Iterable[int],
    kappa: int,
    s_target: int,
    use_codetables: bool = False,
    codetables_timeout: float = 15.0,
) -> tuple[list[CodeCandidate], list[CodeCandidate]]:
    all_rows: list[CodeCandidate] = []
    selected: list[CodeCandidate] = []
    for delta in deltas:
        for r in divisors(delta):
            candidates = all_candidates_for(
                delta,
                kappa,
                s_target,
                r,
                use_codetables=use_codetables,
                codetables_timeout=codetables_timeout,
            )
            all_rows.extend(candidates)
            selected.append(min(candidates, key=candidate_key))
    return all_rows, selected


def write_csv_rows(path: Path, rows: list[CodeCandidate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "delta",
                "kappa",
                "s",
                "r",
                "q",
                "k",
                "s_r",
                "d_star",
                "method",
                "n",
                "d_guaranteed",
                "n*r",
                "griesmer",
                "gap",
                "overshoot",
                "parent_n",
                "parent_k",
                "bch_start",
                "bch_degree",
                "external_lower_bound",
                "external_upper_bound",
                "source_url",
                "verified",
                "note",
            ]
        )
        for c in rows:
            writer.writerow(
                [
                    c.delta,
                    c.kappa,
                    c.s_target,
                    c.r,
                    c.q,
                    c.k,
                    c.s_r,
                    c.d_star,
                    c.method,
                    c.n_r,
                    c.d_actual,
                    c.bit_proxy,
                    c.griesmer_bound,
                    f"{c.griesmer_gap:.6g}",
                    f"{c.distance_overshoot:.6g}",
                    "" if c.parent_n is None else c.parent_n,
                    "" if c.parent_k is None else c.parent_k,
                    "" if c.bch_start is None else c.bch_start,
                    "" if c.bch_degree is None else c.bch_degree,
                    "" if c.external_lower_bound is None else c.external_lower_bound,
                    "" if c.external_upper_bound is None else c.external_upper_bound,
                    c.source_url,
                    c.verified,
                    c.note,
                ]
            )


def write_latex_table(path: Path, selected: list[CodeCandidate]) -> None:
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Selected code candidates after BCH and best-known-table search for $\kappa=s=128$. The selected length is the shortest among juxtaposed simplex, shortened BCH, exact $k\in\{1,2\}$ Griesmer constructions, and CodeTables witnesses where the public table covers the target range~\cite{GrasslCodeTables,Grassl2006SearchCodes}.}",
        r"\label{tab:code-search}",
        r"\begin{tabular}{rrrrL{.24\linewidth}rrrr}",
        r"\toprule",
        r"$\delta$ & $q$ & $k$ & $s_r$ & method & $n$ & $d$ & $g_q$ & gap \\",
        r"\midrule",
    ]
    for c in selected:
        lines.append(
            f"{c.delta} & {c.q} & {c.k} & {c.s_r} & {tex_escape(c.method)} & "
            f"{c.n_r} & {c.d_actual} & {c.griesmer_bound} & {c.griesmer_gap:.2f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def load_selected_from_csv(path: Path) -> list[CodeCandidate]:
    rows: list[CodeCandidate] = []
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(
                CodeCandidate(
                    delta=int(row["delta"]),
                    kappa=int(row["kappa"]),
                    s_target=int(row["s"]),
                    r=int(row["r"]),
                    q=int(row["q"]),
                    k=int(row["k"]),
                    s_r=int(row["s_r"]),
                    d_star=int(row["d_star"]),
                    method=row["method"],
                    n_r=int(row["n"]),
                    d_actual=int(row["d_guaranteed"]),
                    bit_proxy=int(row["n*r"]),
                    griesmer_bound=int(row["griesmer"]),
                    griesmer_gap=float(row["gap"]),
                    distance_overshoot=float(row["overshoot"]),
                    parent_n=int(row["parent_n"]) if row["parent_n"] else None,
                    parent_k=int(row["parent_k"]) if row["parent_k"] else None,
                    bch_start=int(row["bch_start"]) if row["bch_start"] else None,
                    bch_degree=int(row["bch_degree"]) if row["bch_degree"] else None,
                    external_lower_bound=(
                        int(row["external_lower_bound"]) if row.get("external_lower_bound") else None
                    ),
                    external_upper_bound=(
                        int(row["external_upper_bound"]) if row.get("external_upper_bound") else None
                    ),
                    source_url=row.get("source_url", ""),
                    verified=row["verified"],
                    note=row["note"],
                )
            )
    return rows


def ensure_code_search_artifacts(
    deltas: Iterable[int],
    kappa: int,
    s_target: int,
    force: bool = False,
    use_codetables: bool = False,
    codetables_timeout: float = 15.0,
) -> list[CodeCandidate]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    selected_path = DATA_DIR / "code_search_selected.csv"
    all_path = DATA_DIR / "code_search_all.csv"
    table_path = GENERATED_DIR / "code_search_table.tex"
    metadata_path = DATA_DIR / "code_search_metadata.json"
    if selected_path.exists() and table_path.exists() and not force:
        metadata = {}
        if metadata_path.exists():
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not use_codetables or metadata.get("use_codetables"):
            return load_selected_from_csv(selected_path)
    all_rows, selected = selected_candidates(
        deltas,
        kappa,
        s_target,
        use_codetables=use_codetables,
        codetables_timeout=codetables_timeout,
    )
    write_csv_rows(all_path, all_rows)
    write_csv_rows(selected_path, selected)
    write_latex_table(table_path, selected)
    metadata_path.write_text(
        json.dumps(
            {
                "deltas": list(deltas),
                "kappa": kappa,
                "s": s_target,
                "use_codetables": use_codetables,
                "accessed": date.today().isoformat() if use_codetables else "",
                "note": "CodeTables probes are best-known-table witnesses, not local generator-matrix verification.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description="Search q-ary code candidates.")
    parser.add_argument("--deltas", nargs="+", type=int, default=[4, 8, 12])
    parser.add_argument("--kappa", type=int, default=128)
    parser.add_argument("--s", type=int, default=128)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--online-codetables", action="store_true")
    parser.add_argument("--codetables-timeout", type=float, default=5.0)
    args = parser.parse_args()

    selected = ensure_code_search_artifacts(
        args.deltas,
        args.kappa,
        args.s,
        force=args.force,
        use_codetables=args.online_codetables,
        codetables_timeout=args.codetables_timeout,
    )
    for c in selected:
        print(
            f"delta={c.delta} q={c.q} k={c.k} method={c.method} "
            f"n={c.n_r} d>={c.d_actual} griesmer={c.griesmer_bound} gap={c.griesmer_gap:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
