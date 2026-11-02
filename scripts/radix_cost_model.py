"""Finite-batch model for radix-matched q-ary OT extension.

The model uses the following CCG18 q-ary OT-extension parameters:

  N = 2^delta, q = 2^r, r | delta, k = delta / r
  s_r = r * ceil(s / r) so that the F_2-subfield check condition r | s_r holds
  d*_r = max(kappa, s_r)
  L_r = (2^delta - 1) / (2^r - 1)
  D_r = 2^(delta-r)
  t_r = ceil(d*_r / D_r)
  n_r = t_r * L_r

For q=2^r and the subfield consistency check, CCG18 gives:

  encoding bits = (M + 2s_r) n_r r
  check bits    = 2s_r M + 2s_r n_r r + 2s_r k_r r

Since k_r r = delta, the finite-batch extension communication is:

  B_ext,r(M) = (M + 4s_r)n_r r + 2s_r M + 2s_r delta.

Total in-family communication includes n_r binary base OTs of B_base bits each:

  B_total,r(M) = n_r B_base + B_ext,r(M).
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True)
class RadixParams:
    delta: int
    kappa: int
    s_target: int
    r: int
    q: int
    k: int
    s_r: int
    d_star: int
    L_r: int
    D_r: int
    t_r: int
    n_r: int
    d_actual: int
    bit_proxy: int
    griesmer_bound: int
    griesmer_gap: float
    distance_overshoot: float
    distance_aligned: bool


def divisors(n: int) -> list[int]:
    return [i for i in range(1, n + 1) if n % i == 0]


def padded_statistical_security(s_target: int, r: int) -> int:
    return r * math.ceil(s_target / r)


def griesmer_bound(q: int, k: int, d: int) -> int:
    return sum(math.ceil(d / (q**i)) for i in range(k))


def params_for(delta: int, kappa: int, s_target: int, r: int) -> RadixParams:
    if delta <= 0:
        raise ValueError("delta must be positive")
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    if s_target <= 0:
        raise ValueError("s_target must be positive")
    if r <= 0 or delta % r != 0:
        raise ValueError("r must be a positive divisor of delta")

    q = 2**r
    k = delta // r
    s_r = padded_statistical_security(s_target, r)
    d_star = max(kappa, s_r)
    L_r = (2**delta - 1) // (2**r - 1)
    D_r = 2 ** (delta - r)
    t_r = math.ceil(d_star / D_r)
    n_r = t_r * L_r
    d_actual = t_r * D_r
    bit_proxy = n_r * r
    g_bound = griesmer_bound(q, k, d_star)
    return RadixParams(
        delta=delta,
        kappa=kappa,
        s_target=s_target,
        r=r,
        q=q,
        k=k,
        s_r=s_r,
        d_star=d_star,
        L_r=L_r,
        D_r=D_r,
        t_r=t_r,
        n_r=n_r,
        d_actual=d_actual,
        bit_proxy=bit_proxy,
        griesmer_bound=g_bound,
        griesmer_gap=n_r / g_bound,
        distance_overshoot=d_actual / d_star,
        distance_aligned=(d_actual == d_star),
    )


def all_params(deltas: Iterable[int], kappas: Iterable[int], s_values: Iterable[int]) -> list[RadixParams]:
    rows: list[RadixParams] = []
    for delta in deltas:
        for kappa in kappas:
            for s_target in s_values:
                for r in divisors(delta):
                    rows.append(params_for(delta, kappa, s_target, r))
    return rows


def encoding_bits(p: RadixParams, batch_size: float) -> float:
    return (batch_size + 2 * p.s_r) * p.bit_proxy


def check_bits(p: RadixParams, batch_size: float) -> float:
    return 2 * p.s_r * batch_size + 2 * p.s_r * p.bit_proxy + 2 * p.s_r * p.delta


def extension_bits(p: RadixParams, batch_size: float) -> float:
    return encoding_bits(p, batch_size) + check_bits(p, batch_size)


def total_bits(p: RadixParams, batch_size: float, base_ot_bits: float) -> float:
    return p.n_r * base_ot_bits + extension_bits(p, batch_size)


def griesmer_total_lower_bound(p: RadixParams, batch_size: float, base_ot_bits: float) -> float:
    """Communication lower bound for radix p from the Griesmer length bound."""

    g = p.griesmer_bound
    return (
        g * base_ot_bits
        + (batch_size + 4 * p.s_r) * g * p.r
        + 2 * p.s_r * batch_size
        + 2 * p.s_r * p.delta
    )


def best_radix(
    params: Sequence[RadixParams],
    batch_size: float,
    base_ot_bits: float,
) -> RadixParams:
    if not params:
        raise ValueError("params must not be empty")
    return min(
        params,
        key=lambda p: (total_bits(p, batch_size, base_ot_bits), p.n_r, p.r),
    )


def affine_difference_against_binary(
    p: RadixParams,
    binary: RadixParams,
    base_ot_bits: float,
) -> tuple[float, float]:
    """Return intercept and slope for B_p(M)-B_binary(M)."""

    intercept = (
        (p.n_r - binary.n_r) * base_ot_bits
        + 4 * p.s_r * p.bit_proxy
        + 2 * p.s_r * p.delta
        - 4 * binary.s_r * binary.bit_proxy
        - 2 * binary.s_r * binary.delta
    )
    slope = (p.bit_proxy + 2 * p.s_r) - (binary.bit_proxy + 2 * binary.s_r)
    return intercept, slope


def activation_base_ot_threshold(p: RadixParams, binary: RadixParams) -> float | None:
    """Return B_crit where p starts with a lower intercept than binary.

    The value is defined for the common case n_p < n_1. If the returned value is
    negative, then p is already intercept-better for every nonnegative
    B_base.
    """

    denominator = binary.n_r - p.n_r
    if denominator <= 0:
        return None
    numerator = (
        4 * p.s_r * p.bit_proxy
        + 2 * p.s_r * p.delta
        - 4 * binary.s_r * binary.bit_proxy
        - 2 * binary.s_r * binary.delta
    )
    return numerator / denominator


def break_even_against_binary(
    p: RadixParams,
    binary: RadixParams,
    base_ot_bits: float,
) -> tuple[str, float | None]:
    """Return the M region where p beats binary under exact finite-batch bits."""

    intercept, slope = affine_difference_against_binary(p, binary, base_ot_bits)
    if slope == 0:
        return ("all M", None) if intercept < 0 else ("never", None)
    threshold = -intercept / slope
    if slope > 0:
        if threshold <= 0:
            return "never", threshold
        return "M < threshold", threshold
    if intercept <= 0:
        return "all M", threshold
    return "M > threshold", threshold


def integer_region(status: str, threshold: float | None) -> str:
    if status == "all M":
        return "all"
    if status == "never" or threshold is None:
        return "never"
    if status == "M > threshold":
        return f">{math.floor(threshold)}"
    upper = math.ceil(threshold) - 1
    return str(max(0, upper))


def write_params_csv(rows: Sequence[RadixParams], path: Path) -> None:
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
            ]
        )
        for p in rows:
            writer.writerow(
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
                    f"{p.griesmer_gap:.6g}",
                    f"{p.distance_overshoot:.6g}",
                    p.distance_aligned,
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute q-ary OT radix parameters.")
    parser.add_argument("--deltas", nargs="+", type=int, default=[4, 8, 12])
    parser.add_argument("--kappa", nargs="+", type=int, default=[128])
    parser.add_argument("--s", nargs="+", type=int, default=[128])
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    rows = all_params(args.deltas, args.kappa, args.s)
    if args.out:
        write_params_csv(rows, args.out)
    else:
        writer = csv.writer(__import__("sys").stdout)
        writer.writerow(
            ["delta", "kappa", "s", "r", "q", "k", "s_r", "d_star", "n_r", "d_actual", "n_r*r", "griesmer", "gap", "overshoot"]
        )
        for p in rows:
            writer.writerow(
                [
                    p.delta,
                    p.kappa,
                    p.s_target,
                    p.r,
                    p.q,
                    p.k,
                    p.s_r,
                    p.d_star,
                    p.n_r,
                    p.d_actual,
                    p.bit_proxy,
                    p.griesmer_bound,
                    f"{p.griesmer_gap:.6g}",
                    f"{p.distance_overshoot:.6g}",
                ]
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
