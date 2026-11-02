"""Validate representative finite-batch parameter and communication results."""

from __future__ import annotations

from radix_cost_model import break_even_against_binary, integer_region, params_for


EXPECTED_DELTA8 = {
    1: (2, 255, 255),
    2: (4, 170, 340),
    4: (16, 136, 544),
    8: (256, 128, 1024),
}


def main() -> int:
    delta = 8
    kappa = 128
    s = 128
    for r, (q, expected_n, expected_bits) in EXPECTED_DELTA8.items():
        p = params_for(delta, kappa, s, r)
        assert p.q == q, (r, p.q, q)
        assert p.s_r == s, (q, p.s_r, s)
        assert p.n_r == expected_n, (q, p.n_r, expected_n)
        assert p.bit_proxy == expected_bits, (q, p.bit_proxy, expected_bits)
        print(
            f"OK delta={delta} kappa={kappa} s={s} q={p.q}: "
            f"n_r={p.n_r}, n_r*r={p.bit_proxy}, d_actual={p.d_actual}"
        )

    binary = params_for(8, 128, 128, 1)
    q4 = params_for(8, 128, 128, 2)
    status, threshold = break_even_against_binary(q4, binary, 100)
    assert status == "never", (status, threshold)
    print(
        "OK finite-batch comparison: "
        f"delta=8 q=4 base_ot_bits=100 region={integer_region(status, threshold)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
