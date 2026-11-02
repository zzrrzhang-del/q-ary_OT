"""Generate paper tables and figures for the finite-batch model."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

from radix_cost_model import (
    RadixParams,
    activation_base_ot_threshold,
    all_params,
    best_radix,
    break_even_against_binary,
    divisors,
    extension_bits,
    griesmer_total_lower_bound,
    griesmer_bound,
    integer_region,
    params_for,
    padded_statistical_security,
    total_bits,
)
from code_search import CodeCandidate, ensure_code_search_artifacts


ROOT = Path(__file__).resolve().parents[1]
TEX_DIR = ROOT / "tex"
GENERATED_DIR = TEX_DIR / "generated"
FIGURE_DIR = TEX_DIR / "figures"
DATA_DIR = GENERATED_DIR / "data"

DELTAS = [4, 8, 12]
KAPPAS = [128]
S_VALUES = [80, 128, 192, 256]
MAIN_KAPPA = 128
MAIN_S = 128
BASE_OT_BITS = [10**2, 10**3, 10**4, 10**5, 10**6]
MAIN_DELTA = 8
PADDING_S_VALUES = [40, 64, 128]
PADDING_POINTS = [
    (8, 10**3, 100),
    (8, 10**3, 10_000),
    (8, 10**5, 100),
    (8, 10**5, 10_000),
    (12, 10**4, 100),
    (12, 10**4, 10_000),
    (12, 10**5, 100),
    (12, 10**5, 10_000),
]

APPLICATION_CASES = [
    ("4-bit LUT backend", 4, 10**2),
    ("8-bit LUT interactive backend", 8, 10**2),
    ("8-bit LUT moderate backend", 8, 10**4),
    ("12-bit LUT moderate backend", 12, 10**4),
    ("Private sampling backend", 8, 10**2),
    ("Large-batch backend", 8, 10**6),
]


def selected_code_params(force: bool = False) -> list[CodeCandidate]:
    return ensure_code_search_artifacts(
        DELTAS,
        MAIN_KAPPA,
        MAIN_S,
        force=force,
        use_codetables=True,
    )


def selected_for_delta(delta: int, selected: list[CodeCandidate]) -> list[CodeCandidate]:
    return sorted([c for c in selected if c.delta == delta], key=lambda c: c.r)


def ensure_dirs() -> None:
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def tex_escape(value: object) -> str:
    return (
        str(value)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("&", r"\&")
        .replace("%", r"\%")
    )


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def format_region(status: str, threshold: float | None) -> str:
    return integer_region(status, threshold)


def format_number(value: float | None) -> str:
    if value is None:
        return "-"
    if value < 0:
        return "<0"
    if value >= 100000:
        return f"{value:.3g}"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 1:
        return f"{value:.2f}"
    return f"{value:.3g}"


def format_certificate_cell(label: str, cert_q: int | None, selected_q: int) -> str:
    if label == "certified non-binary advantage":
        return f"cert. non-bin. q={cert_q}"
    if label == "certified binary":
        return "cert. binary"
    return f"dep. (selected q={selected_q})"


def compact_certificate_tex(value: str) -> str:
    if value == "cert. binary":
        return "B"
    if value.startswith("cert. non-bin. q="):
        q = value.rsplit("=", 1)[1]
        return rf"NB ($q={q}$)"
    q = value.rsplit("=", 1)[1].rstrip(")")
    return rf"CD ($q_{{\rm sel}}={q}$)"


def critical_batch_size(p: CodeCandidate, binary: CodeCandidate, base_ot_bits: float) -> float | None:
    status, threshold = break_even_against_binary(p, binary, base_ot_bits)
    if status == "M < threshold" and threshold is not None and threshold > 0:
        return threshold
    return None


def certified_region(params: list[CodeCandidate], batch_size: float, base_ot_bits: float) -> tuple[str, int | None, int]:
    binary = next(p for p in params if p.r == 1)
    binary_witness = total_bits(binary, batch_size, base_ot_bits)
    binary_lower = griesmer_total_lower_bound(binary, batch_size, base_ot_bits)
    nonbinary_certified = [
        p
        for p in params
        if p.r > 1 and total_bits(p, batch_size, base_ot_bits) < binary_lower
    ]
    if nonbinary_certified:
        best = min(nonbinary_certified, key=lambda p: (total_bits(p, batch_size, base_ot_bits), p.r))
        return "certified non-binary advantage", best.q, best_radix(params, batch_size, base_ot_bits).q
    if all(
        griesmer_total_lower_bound(p, batch_size, base_ot_bits) >= binary_witness
        for p in params
        if p.r > 1
    ):
        return "certified binary", 2, best_radix(params, batch_size, base_ot_bits).q
    return "code-dependent", None, best_radix(params, batch_size, base_ot_bits).q


def pareto_frontier(params: list[RadixParams]) -> list[RadixParams]:
    frontier: list[RadixParams] = []
    for p in sorted(params, key=lambda x: (x.n_r, x.bit_proxy)):
        dominated = False
        for other in params:
            if other == p:
                continue
            if (
                other.n_r <= p.n_r
                and other.bit_proxy <= p.bit_proxy
                and (other.n_r < p.n_r or other.bit_proxy < p.bit_proxy)
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(p)
    return frontier


def adjusted_candidate_for_s(candidate: CodeCandidate, s_target: int) -> CodeCandidate:
    """Reuse selected witnesses while recomputing padding-sensitive quantities."""

    s_r = padded_statistical_security(s_target, candidate.r)
    d_star = max(candidate.kappa, s_r)
    g_bound = griesmer_bound(candidate.q, candidate.k, d_star)
    n_r = candidate.n_r
    d_actual = candidate.d_actual
    method = candidate.method
    verified = candidate.verified
    note = candidate.note

    if candidate.k == 1:
        n_r = d_star
        d_actual = d_star
        method = "repetition"
        verified = "explicit"
        note = "optimal for k=1; recomputed for padding-sensitivity sweep"
    elif candidate.k == 2 and candidate.method in {"projective-line Griesmer", "repetition"}:
        n_r = g_bound
        d_actual = d_star
        method = "projective-line Griesmer"
        verified = "theorem"
        note = "Griesmer-optimal k=2 witness recomputed for padding-sensitivity sweep"

    if d_actual < d_star:
        raise ValueError(
            f"Selected witness q={candidate.q} has d={d_actual}, below adjusted target {d_star}"
        )

    return replace(
        candidate,
        s_target=s_target,
        s_r=s_r,
        d_star=d_star,
        method=method,
        n_r=n_r,
        d_actual=d_actual,
        bit_proxy=n_r * candidate.r,
        griesmer_bound=g_bound,
        griesmer_gap=n_r / g_bound,
        distance_overshoot=d_actual / d_star,
        verified=verified,
        note=note,
    )


def selected_for_s(selected: list[CodeCandidate], s_target: int) -> list[CodeCandidate]:
    return [adjusted_candidate_for_s(candidate, s_target) for candidate in selected]


def params_rows() -> list[RadixParams]:
    return all_params(DELTAS, KAPPAS, S_VALUES)


def write_params_artifacts() -> None:
    rows = params_rows()
    write_csv(
        DATA_DIR / "params_all.csv",
        [
            "delta",
            "kappa",
            "s",
            "r",
            "q",
            "k",
            "s_r",
            "d_star",
            "L_r",
            "D_r",
            "t_r",
            "n_r",
            "d_actual",
            "n_r*r",
            "griesmer",
            "gap",
            "overshoot",
            "aligned",
        ],
        [
            [
                p.delta,
                p.kappa,
                p.s_target,
                p.r,
                p.q,
                p.k,
                p.s_r,
                p.d_star,
                p.L_r,
                p.D_r,
                p.t_r,
                p.n_r,
                p.d_actual,
                p.bit_proxy,
                p.griesmer_bound,
                f"{p.griesmer_gap:.3f}",
                f"{p.distance_overshoot:.3f}",
                p.distance_aligned,
            ]
            for p in rows
        ],
    )

    main_rows = [
        params_for(delta, MAIN_KAPPA, MAIN_S, r)
        for delta in DELTAS
        for r in divisors(delta)
    ]
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Repeated-simplex construction parameters for $\kappa=s=128$. The columns $g_q$ and gap report the Griesmer lower bound and the ratio $n_r/g_q(k,d^\star_r)$.}",
        r"\label{tab:params}",
        r"\begin{tabular}{rrrrrrrrrr}",
        r"\toprule",
        r"$\delta$ & $q$ & $s_r$ & $d^\star_r$ & $n_r$ & $d_{\rm act}$ & $n_r r$ & $g_q$ & gap & over. \\",
        r"\midrule",
    ]
    for p in main_rows:
        lines.append(
            f"{p.delta} & {p.q} & {p.s_r} & {p.d_star} & {p.n_r} & "
            f"{p.d_actual} & {p.bit_proxy} & {p.griesmer_bound} & "
            f"{p.griesmer_gap:.2f} & {p.distance_overshoot:.2f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    (GENERATED_DIR / "params_table.tex").write_text("\n".join(lines), encoding="utf-8")


def verification_mark(candidate: CodeCandidate) -> str:
    if candidate.delta == 12 and candidate.q in {2, 8, 16}:
        return "Magma"
    if candidate.delta == 12 and candidate.q == 4:
        return "CT+Magma"
    if candidate.method == "juxtaposed simplex":
        return "formula"
    if candidate.method == "projective-line Griesmer":
        return "theorem"
    if candidate.method == "repetition":
        return "explicit"
    return candidate.verified or "witness"


def short_method_name(method: str) -> str:
    if method == "juxtaposed simplex":
        return "simplex"
    if method == "projective-line Griesmer":
        return "proj.-line"
    if method == "shortened BCH":
        return "BCH"
    if method == "CodeTables witness":
        return "CodeTables"
    return method


def write_construction_witness_artifacts(selected: list[CodeCandidate]) -> None:
    rows: list[list[object]] = []
    for candidate in sorted(selected, key=lambda c: (c.delta, c.r)):
        simplex = params_for(candidate.delta, MAIN_KAPPA, MAIN_S, candidate.r)
        rows.append(
            [
                candidate.delta,
                candidate.q,
                candidate.k,
                candidate.s_r,
                candidate.d_star,
                simplex.n_r,
                candidate.n_r,
                candidate.griesmer_bound,
                candidate.method,
                verification_mark(candidate),
            ]
        )
    write_csv(
        DATA_DIR / "construction_witnesses.csv",
        [
            "delta",
            "q",
            "k",
            "s_r",
            "d_star",
            "simplex_n",
            "selected_n",
            "griesmer",
            "method",
            "verification",
        ],
        rows,
    )

    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Construction witnesses for $\kappa=s=128$. The column $n_{\rm simp}$ is the repeated-simplex length, while $n$ is the selected admissible witness used in the cost model. The check column records the proof or verification source.}",
        r"\label{tab:construction-witnesses}",
        r"\begingroup",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular}{rrrrrrrll}",
        r"\toprule",
        r"$\delta$ & $q$ & $k$ & $s_r$ & $d^\star$ & $n_{\rm simp}$ & $n$ & method & check \\",
        r"\midrule",
    ]
    for delta, q, k, s_r, d_star, simplex_n, selected_n, _griesmer, method, check in rows:
        lines.append(
            f"{delta} & {q} & {k} & {s_r} & {d_star} & {simplex_n} & {selected_n} & "
            f"{tex_escape(short_method_name(str(method)))} & {tex_escape(check)} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\endgroup", r"\end{table}", ""])
    (GENERATED_DIR / "construction_witness_table.tex").write_text("\n".join(lines), encoding="utf-8")


def write_symbolic_points_table() -> None:
    rows = [
        [
            "sub-classical stress",
            10**2,
            "",
            "below reported public-key base-OT communication; included only as a lower stress point",
        ],
        [
            "classical EOT anchor",
            10**3,
            "CO15 112 B, Endemic-DDH 145 B, NP01 165 B",
            "same order as reported classical one-out-of-two base OTs",
        ],
        [
            "intermediate",
            10**4,
            "",
            "between classical and Kyber-sized base-OT messages",
        ],
        [
            "Kyber-scale sensitivity point",
            10**5,
            "Endemic-Kyber 5,934 B",
            "conservative point above the reported Kyber-based base OT",
        ],
        [
            "stress upper point",
            10**6,
            "",
            "deployment stress point for extra framing, latency, or constrained-client overhead",
        ],
    ]
    write_csv(
        DATA_DIR / "base_ot_points.csv",
        ["label", "B_base_bits", "measurement_anchor", "interpretation"],
        rows,
    )
    stale_table = GENERATED_DIR / "base_ot_points_table.tex"
    if stale_table.exists():
        stale_table.unlink()


def write_breakeven_artifacts() -> None:
    selected = selected_code_params()
    rows: list[list[object]] = []
    for delta in DELTAS:
        params = selected_for_delta(delta, selected)
        binary = next(p for p in params if p.r == 1)
        for base_bits in BASE_OT_BITS:
            for p in params:
                if p.r == 1:
                    continue
                status, threshold = break_even_against_binary(p, binary, base_bits)
                rows.append(
                    [
                        delta,
                        MAIN_KAPPA,
                        MAIN_S,
                        base_bits,
                        p.q,
                        status,
                        "" if threshold is None else f"{threshold:.6g}",
                        format_region(status, threshold),
                    ]
                )
    write_csv(
        DATA_DIR / "breakeven_all.csv",
        ["delta", "kappa", "s", "B_base_bits", "q", "status", "threshold", "integer_region"],
        rows,
    )

    delta8_rows = [row for row in rows if row[0] == MAIN_DELTA]
    by_bits: dict[int, dict[int, str]] = {bits: {} for bits in BASE_OT_BITS}
    for _delta, _kappa, _s, bits, q, _status, _threshold, region in delta8_rows:
        by_bits[int(bits)][int(q)] = str(region)
    q_values = [4, 16, 256]
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\caption{Largest integer batch size $M$ for which $q>2$ beats the in-family binary baseline for $\delta=8$, $\kappa=s=128$, under exact finite-batch communication.}",
        r"\label{tab:breakeven}",
        r"\begin{tabular}{rrrr}",
        r"\toprule",
        r"$B_{\rm base}$ bits & $q=4$ & $q=16$ & $q=256$ \\",
        r"\midrule",
    ]
    for bits in BASE_OT_BITS:
        cells = [by_bits[bits].get(q, "-") for q in q_values]
        lines.append(f"{bits} & {cells[0]} & {cells[1]} & {cells[2]} \\\\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    breakeven_table = GENERATED_DIR / "breakeven_table.tex"
    breakeven_table.write_text("\n".join(lines), encoding="utf-8")
    breakeven_table.unlink()

    summary_rows: list[list[object]] = []
    for delta in DELTAS:
        params = selected_for_delta(delta, selected)
        binary = next(p for p in params if p.r == 1)
        for bits in BASE_OT_BITS:
            beating = []
            for p in params:
                if p.r == 1:
                    continue
                status, threshold = break_even_against_binary(p, binary, bits)
                if status != "never":
                    beating.append(f"q={p.q}:{format_region(status, threshold)}")
            best_at_100 = best_radix(params, 100, bits)
            best_at_10000 = best_radix(params, 10000, bits)
            summary_rows.append(
                [
                    delta,
                    bits,
                    "; ".join(beating) if beating else "none",
                    best_at_100.q,
                    best_at_10000.q,
                ]
            )
    write_csv(
        DATA_DIR / "breakeven_summary.csv",
        ["delta", "B_base_bits", "radices_beating_binary", "best_q_M100", "best_q_M10000"],
        summary_rows,
    )
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Summary of non-binary radices beating the in-family binary baseline. Entries list exact finite-batch regions; $q^\star$ columns show the global minimum-cost radix at two representative batch sizes.}",
        r"\label{tab:breakeven-summary}",
        r"\begin{tabular}{rrL{.42\linewidth}rr}",
        r"\toprule",
        r"$\delta$ & $B_{\rm base}$ & beating regions & $q^\star_{M=100}$ & $q^\star_{M=10^4}$ \\",
        r"\midrule",
    ]
    for delta, bits, beating, best100, best10000 in summary_rows:
        lines.append(f"{delta} & {bits} & {tex_escape(beating)} & {best100} & {best10000} \\\\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    breakeven_summary_table = GENERATED_DIR / "breakeven_summary_table.tex"
    breakeven_summary_table.write_text("\n".join(lines), encoding="utf-8")
    breakeven_summary_table.unlink()


def write_threshold_artifacts(selected: list[CodeCandidate]) -> None:
    threshold_rows: list[list[object]] = []
    for delta in DELTAS:
        params = selected_for_delta(delta, selected)
        binary = next(p for p in params if p.r == 1)
        for p in params:
            if p.r == 1:
                continue
            bcrit = activation_base_ot_threshold(p, binary)
            threshold_rows.append(
                [
                    delta,
                    p.q,
                    p.n_r,
                    p.griesmer_bound,
                    "" if bcrit is None else f"{bcrit:.6g}",
                    p.bit_proxy + 2 * p.s_r - (binary.bit_proxy + 2 * binary.s_r),
                    "" if critical_batch_size(p, binary, 10**3) is None else f"{critical_batch_size(p, binary, 10**3):.6g}",
                    "" if critical_batch_size(p, binary, 10**5) is None else f"{critical_batch_size(p, binary, 10**5):.6g}",
                ]
            )
    write_csv(
        DATA_DIR / "critical_thresholds.csv",
        [
            "delta",
            "q",
            "n_witness",
            "griesmer",
            "B_crit_bits",
            "delta_slope",
            "Mcrit_B1000",
            "Mcrit_B100000",
        ],
        threshold_rows,
    )

    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Activation thresholds against the in-family binary witness. $B_{\rm crit}$ is the base-OT communication threshold; $M_{\rm crit}$ is the finite-batch upper threshold after activation. A dash indicates that $B_{\rm base}\leq B_{\rm crit}$, so the non-binary witness is not activated and no positive $M_{\rm crit}$ exists.}",
        r"\label{tab:critical-thresholds}",
        r"\begin{tabular}{rrrrr}",
        r"\toprule",
        r"$\delta$ & $q$ & $B_{\rm crit}$ & $M_{\rm crit}(10^3)$ & $M_{\rm crit}(10^5)$ \\",
        r"\midrule",
    ]
    for delta, q, n, g, bcrit, _slope, m1000, m100000 in threshold_rows:
        lines.append(
            f"{delta} & {q} & {format_number(float(bcrit)) if bcrit else '-'} & "
            f"{format_number(float(m1000)) if m1000 else '-'} & {format_number(float(m100000)) if m100000 else '-'} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    (GENERATED_DIR / "critical_thresholds_table.tex").write_text("\n".join(lines), encoding="utf-8")


def write_certification_point_artifacts(selected: list[CodeCandidate]) -> None:
    rows: list[list[object]] = []
    for delta in DELTAS:
        params = selected_for_delta(delta, selected)
        for bits in BASE_OT_BITS:
            row: list[object] = [delta, bits]
            for batch in [100, 10_000]:
                label, cert_q, selected_q = certified_region(params, batch, bits)
                if label == "certified non-binary advantage":
                    cell = f"certified non-binary q={cert_q}"
                elif label == "certified binary":
                    cell = "certified binary"
                else:
                    cell = f"code-dependent (selected q={selected_q})"
                row.append(cell)
            rows.append(row)
    write_csv(
        DATA_DIR / "certification_points.csv",
        ["delta", "B_base_bits", "M100_class", "M10000_class"],
        rows,
    )
    stale_table = GENERATED_DIR / "certification_points_table.tex"
    if stale_table.exists():
        stale_table.unlink()


def write_padding_sensitivity_artifacts(selected: list[CodeCandidate]) -> None:
    rows: list[list[object]] = []
    for s_target in PADDING_S_VALUES:
        adjusted = selected_for_s(selected, s_target)
        for delta, bits, batch in PADDING_POINTS:
            params = selected_for_delta(delta, adjusted)
            label, cert_q, selected_q = certified_region(params, batch, bits)
            rows.append(
                [
                    s_target,
                    delta,
                    bits,
                    batch,
                    format_certificate_cell(label, cert_q, selected_q),
                    selected_q,
                ]
            )
    write_csv(
        DATA_DIR / "padding_sensitivity.csv",
        ["s", "delta", "B_base_bits", "M", "certificate", "selected_q"],
        rows,
    )

    row_keys = PADDING_POINTS
    by_key = {
        (int(s), int(delta), int(bits), int(batch)): str(cert)
        for s, delta, bits, batch, cert, _selected_q in rows
    }
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Padding-sensitivity check at $\kappa=128$. B, NB, and CD denote certified binary, certified non-binary advantage, and code-dependent, respectively. For the available constructions, the sampled certificate class is stable across $s\in\{40,64,128\}$; one high-cost $\delta=12$ point changes the selected non-binary construction.}",
        r"\label{tab:padding-sensitivity}",
        r"\begingroup",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular}{rrrL{.19\linewidth}L{.19\linewidth}L{.19\linewidth}}",
        r"\toprule",
        r"$\delta$ & $B_{\rm base}$ & $M$ & $s=40$ & $s=64$ & $s=128$ \\",
        r"\midrule",
    ]
    for delta, bits, batch in row_keys:
        cells = [
            by_key[(s_target, delta, bits, batch)]
            for s_target in PADDING_S_VALUES
        ]
        lines.append(
            f"{delta} & {bits} & {batch} & "
            f"{compact_certificate_tex(cells[0])} & {compact_certificate_tex(cells[1])} & {compact_certificate_tex(cells[2])} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\endgroup", r"\end{table}", ""])
    (GENERATED_DIR / "padding_sensitivity_table.tex").write_text("\n".join(lines), encoding="utf-8")


def write_tradeoff_plot(delta: int, selected: list[CodeCandidate]) -> None:
    params = selected_for_delta(delta, selected)
    frontier = pareto_frontier(params)
    write_csv(
        DATA_DIR / f"pareto_delta{delta}.csv",
        ["q", "r", "n_r", "n_r*r", "pareto", "griesmer_gap", "overshoot"],
        [
            [
                p.q,
                p.r,
                p.n_r,
                p.bit_proxy,
                p in frontier,
                f"{p.griesmer_gap:.3f}",
                f"{p.distance_overshoot:.3f}",
            ]
            for p in params
        ],
    )

    fig, ax = plt.subplots(figsize=(5.3, 3.6))
    for p in params:
        is_front = p in frontier
        ax.scatter(
            p.n_r,
            p.bit_proxy,
            s=75 if is_front else 45,
            facecolors="tab:blue" if is_front else "none",
            edgecolors="tab:blue" if is_front else "tab:gray",
        )
        ax.annotate(f"q={p.q}", (p.n_r, p.bit_proxy), xytext=(5, 5), textcoords="offset points", fontsize=8)
    frontier_sorted = sorted(frontier, key=lambda p: p.n_r)
    ax.plot([p.n_r for p in frontier_sorted], [p.bit_proxy for p in frontier_sorted], "--", linewidth=1)
    ax.set_xlabel("Base OT count $n_r$")
    ax.set_ylabel("Encoding bit coefficient $n_r r$")
    ax.set_title(rf"Selected-code trade-off for $\delta={delta}$")
    ax.grid(True, linewidth=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"pareto_delta{delta}.pdf")
    plt.close(fig)


def write_tradeoff_overview_plot(selected: list[CodeCandidate]) -> None:
    fig, axes = plt.subplots(1, len(DELTAS), figsize=(7.4, 3.35), sharey=False)
    for ax, delta in zip(axes, DELTAS):
        params = selected_for_delta(delta, selected)
        frontier = pareto_frontier(params)
        for p in params:
            is_front = p in frontier
            ax.scatter(
                p.n_r,
                p.bit_proxy,
                s=65 if is_front else 38,
                facecolors="#2b6cb0" if is_front else "white",
                edgecolors="#2b6cb0" if is_front else "#777777",
                linewidths=1.1,
            )
            ax.annotate(f"{p.q}", (p.n_r, p.bit_proxy), xytext=(4, 4), textcoords="offset points", fontsize=7)
        frontier_sorted = sorted(frontier, key=lambda p: p.n_r)
        ax.plot([p.n_r for p in frontier_sorted], [p.bit_proxy for p in frontier_sorted], "--", linewidth=0.9, color="#555555")
        ax.set_title(rf"$\delta={delta}$", fontsize=9)
        ax.set_xlabel("$n_r$", fontsize=8)
        ax.grid(True, linewidth=0.35, alpha=0.45)
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("$n_r r$", fontsize=8)
    fig.suptitle("Selected-code Pareto view", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIGURE_DIR / "pareto_overview.pdf", bbox_inches="tight")
    plt.close(fig)


def write_region_plot(delta: int, selected: list[CodeCandidate]) -> None:
    params = selected_for_delta(delta, selected)
    bits_grid = np.logspace(2, 6, 181)
    m_grid = np.logspace(0, 6, 181)
    z = np.zeros((len(m_grid), len(bits_grid)), dtype=int)
    labels = ["certified binary", "code-dependent", "certified non-binary advantage"]
    label_to_index = {label: idx for idx, label in enumerate(labels)}
    rows: list[list[object]] = []
    for i, m in enumerate(m_grid):
        for j, bits in enumerate(bits_grid):
            label, cert_q, selected_q = certified_region(params, m, bits)
            z[i, j] = label_to_index[label]
        if i in {0, 45, 90, 135, 180}:
            for bits in BASE_OT_BITS:
                label, cert_q, selected_q = certified_region(params, m, bits)
                rows.append([f"{m:.6g}", bits, label, "" if cert_q is None else cert_q, selected_q])
    write_csv(
        DATA_DIR / f"radix_regions_delta{delta}_sample.csv",
        ["M", "B_base_bits", "certificate_class", "certified_nonbinary_q", "selected_q"],
        rows,
    )

    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    cmap = ListedColormap(["#4c78a8", "#f2c14e", "#59a14f"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    image = ax.imshow(
        z,
        origin="lower",
        extent=(2, 6, 0, 6),
        interpolation="nearest",
        aspect="auto",
        cmap=cmap,
        norm=norm,
    )
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_ticks(range(len(labels)))
    cbar.set_ticklabels(labels)
    cbar.set_label("certificate class")
    ax.set_xticks([2, 3, 4, 5, 6])
    ax.set_xticklabels([r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"])
    ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    ax.set_yticklabels([r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"])
    ax.set_xlabel(r"$B_{\rm base}$ bits")
    ax.set_ylabel("Batch size M")
    ax.set_title(rf"Certified design regions for $\delta={delta}$")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"radix_regions_delta{delta}.pdf")
    plt.close(fig)


def write_region_overview_plot(selected: list[CodeCandidate]) -> None:
    bits_grid = np.logspace(2, 6, 181)
    m_grid = np.logspace(0, 6, 181)
    labels = ["certified binary", "code-dependent", "certified non-binary advantage"]
    label_to_index = {label: idx for idx, label in enumerate(labels)}
    black = "#000000"
    white = "#ffffff"
    blue = "#1f5a94"
    red = "#b22222"
    category_colors = [black, blue, red]
    category_hatches = ["", "////", "\\\\\\\\"]
    x_grid = np.linspace(2, 6, len(bits_grid))
    y_grid = np.linspace(0, 6, len(m_grid))

    fig, axes = plt.subplots(len(DELTAS), 1, figsize=(6.5, 7.6), sharex=False, sharey=True)
    for ax, delta in zip(axes, DELTAS):
        params = selected_for_delta(delta, selected)
        z = np.zeros((len(m_grid), len(bits_grid)), dtype=int)
        for i, batch_size in enumerate(m_grid):
            for j, base_bits in enumerate(bits_grid):
                label, _cert_q, _selected_q = certified_region(params, batch_size, base_bits)
                z[i, j] = label_to_index[label]
        for category_index, (color, hatch) in enumerate(zip(category_colors, category_hatches)):
            if not np.any(z == category_index):
                continue
            masked = np.ma.masked_where(z != category_index, z)
            region = ax.contourf(
                x_grid,
                y_grid,
                masked,
                levels=[category_index - 0.5, category_index + 0.5],
                colors=[white],
                hatches=[hatch],
                antialiased=True,
            )
            region.set_edgecolor(color)
            region.set_linewidth(0.25)
            boundary = (z == category_index).astype(float)
            if np.any(boundary == 0):
                ax.contour(x_grid, y_grid, boundary, levels=[0.5], colors=black, linewidths=0.7)
        ax.set_ylabel(rf"$M$ ($\delta={delta}$)", fontsize=8.5)
        ax.set_xlabel(r"$B_{\rm base}$ bits", fontsize=8.5)
        ax.set_xticks([2, 3, 4, 5, 6])
        ax.set_xticklabels([r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"])
        ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
        ax.set_yticklabels([r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"])
        ax.grid(False)
        ax.tick_params(labelsize=7.5, colors=black)
        for spine in ax.spines.values():
            spine.set_color(black)
    legend_handles = [
        Patch(facecolor=white, edgecolor=black, label="certified binary"),
        Patch(facecolor=white, edgecolor=blue, hatch=category_hatches[1], label="code-dependent"),
        Patch(facecolor=white, edgecolor=red, hatch=category_hatches[2], label="certified non-binary advantage"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.50, 0.015),
        frameon=False,
        fontsize=7.2,
        handlelength=2.2,
        ncol=3,
        columnspacing=1.2,
    )
    fig.suptitle("Certified CCG-family design regions", fontsize=10, color=black)
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.14, top=0.91, hspace=0.42)
    fig.savefig(FIGURE_DIR / "radix_regions_overview.pdf", bbox_inches="tight")
    plt.close(fig)


def write_analysis_framework_figure() -> None:
    black = "#000000"
    white = "#ffffff"
    fig, ax = plt.subplots(figsize=(7.2, 2.75))
    ax.set_axis_off()
    blocks = [
        (0.02, 0.56, 0.21, 0.27, "Public workload\n$\\delta, M, \\kappa, s, B_{base}$", white, black, black),
        (0.27, 0.56, 0.21, 0.27, "Admissible CCG\nradix/code family", white, black, black),
        (0.52, 0.56, 0.21, 0.27, "Exact affine cost\n+ Griesmer bounds", white, black, black),
        (0.77, 0.56, 0.21, 0.27, "Certified\nworkload regions", white, black, black),
        (0.11, 0.13, 0.23, 0.23, "Code witnesses\nsimplex / BCH / tables", white, black, black),
        (0.39, 0.13, 0.23, 0.23, "Grounding checks\nbase-OT anchors\npadding sweep", white, black, black),
        (0.67, 0.13, 0.23, 0.23, "Main result\nFig. 3 phase diagram", white, black, black),
    ]
    for x, y, w, h, text, face_color, edge_color, text_color in blocks:
        ax.add_patch(Rectangle((x, y), w, h, linewidth=1.1, edgecolor=edge_color, facecolor=face_color))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.3, color=text_color)
    arrows = [
        ((0.23, 0.695), (0.27, 0.695), black),
        ((0.48, 0.695), (0.52, 0.695), black),
        ((0.73, 0.695), (0.77, 0.695), black),
        ((0.34, 0.245), (0.39, 0.245), black),
        ((0.62, 0.245), (0.67, 0.245), black),
        ((0.215, 0.36), (0.35, 0.56), black),
        ((0.505, 0.36), (0.625, 0.56), black),
        ((0.875, 0.56), (0.785, 0.36), black),
    ]
    for start, end, color in arrows:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10, linewidth=1.0, color=color))
    ax.text(
        0.50,
        0.93,
        "Analysis framework: from public workload to certified regions",
        ha="center",
        fontsize=10,
        weight="bold",
        color=black,
    )
    fig.savefig(FIGURE_DIR / "analysis_framework.pdf", bbox_inches="tight")
    plt.close(fig)


def write_base_ot_sensitivity_plot() -> None:
    params = selected_for_delta(MAIN_DELTA, selected_code_params())
    binary = next(p for p in params if p.r == 1)
    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    rows: list[list[object]] = []
    for p in params:
        if p.r == 1:
            continue
        values = []
        for bits in BASE_OT_BITS:
            status, threshold = break_even_against_binary(p, binary, bits)
            if status == "M < threshold" and threshold is not None and threshold > 0:
                value = threshold
            elif status == "all M":
                value = 1e6
            else:
                value = 0.0
            values.append(value)
            rows.append([p.q, bits, status, "" if threshold is None else f"{threshold:.6g}", format_region(status, threshold)])
        ax.plot(BASE_OT_BITS, values, marker="o", label=f"q={p.q}")
    write_csv(
        DATA_DIR / "base_ot_sensitivity.csv",
        ["q", "B_base_bits", "status", "threshold", "integer_region"],
        rows,
    )
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel(r"$B_{\rm base}$ bits")
    ax.set_ylabel("Break-even M")
    ax.set_title(r"Exact break-even sensitivity for $\delta=8$")
    ax.grid(True, linewidth=0.4, alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "base_ot_sensitivity.pdf")
    plt.close(fig)


def write_application_artifacts() -> None:
    selected = selected_code_params()
    rows: list[list[object]] = []
    for name, delta, batch in APPLICATION_CASES:
        params = selected_for_delta(delta, selected)
        binary = next(p for p in params if p.r == 1)
        row: list[object] = [name, delta, 2**delta, batch]
        for bits in [10**3, 10**5, 10**6]:
            best = best_radix(params, batch, bits)
            ratio = total_bits(best, batch, bits) / total_bits(binary, batch, bits)
            row.extend([best.q, f"{ratio:.3f}"])
        rows.append(row)
    write_csv(
        DATA_DIR / "application_cases.csv",
        [
            "workload",
            "delta",
            "N",
            "M",
            "B1000_q",
            "B1000_ratio",
            "B100000_q",
            "B100000_ratio",
            "B1000000_q",
            "B1000000_ratio",
        ],
        rows,
    )
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Radix-dependent preprocessing ratios. Ratios compare exact backend communication against the in-family binary baseline and exclude payload masking, online LUT communication, and other application phases. The last column reports the minimum communication ratio to the corresponding binary CCG instance over the displayed $B_{\rm base}$ points.}",
        r"\label{tab:application-cases}",
        r"\begin{tabular}{L{.30\linewidth}rrrrrrr}",
        r"\toprule",
        r"Backend workload & $\delta$ & $N$ & $M$ & $q^\star_{10^3}$ & $q^\star_{10^5}$ & $q^\star_{10^6}$ & \shortstack{min. comm.\\ratio} \\",
        r"\midrule",
    ]
    for row in rows:
        name, delta, n_value, batch = row[:4]
        q_values = [row[4], row[6], row[8]]
        ratios = [float(row[5]), float(row[7]), float(row[9])]
        lines.append(
            f"{tex_escape(name)} & {delta} & {n_value} & {batch} & "
            f"{q_values[0]} & {q_values[1]} & {q_values[2]} & {min(ratios):.3f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    (GENERATED_DIR / "application_cases_table.tex").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ensure_dirs()
    selected = selected_code_params()
    write_analysis_framework_figure()
    write_params_artifacts()
    write_construction_witness_artifacts(selected)
    write_symbolic_points_table()
    write_breakeven_artifacts()
    write_threshold_artifacts(selected)
    write_certification_point_artifacts(selected)
    write_padding_sensitivity_artifacts(selected)
    for delta in DELTAS:
        write_tradeoff_plot(delta, selected)
        write_region_plot(delta, selected)
    write_tradeoff_overview_plot(selected)
    write_region_overview_plot(selected)
    write_base_ot_sensitivity_plot()
    write_application_artifacts()
    print(f"Generated artifacts under {GENERATED_DIR} and {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
