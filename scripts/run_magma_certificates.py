"""Run Magma certificates for selected q-ary code candidates.

This script keeps source code under `scripts/` and writes reproducibility
artifacts under `tex/generated/`. It verifies the construction-level code
parameters used by the paper with Magma whenever the local Magma installation
can construct the code.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import ssl
import subprocess
import tempfile
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED_DIR = ROOT / "tex" / "generated"
DATA_DIR = GENERATED_DIR / "data"
DEFAULT_MAGMA = Path("magma")
CODETABLES_Q4_URL = "https://www.codetables.de/BKLC/BKLC.php?k=6&n=177&q=4"
CODETABLES_Q4_MATRIX_SHA256 = "12252741684441389e88b10aa8c0dfd34858ae9916a6b105aebba349d91727cf"
METHOD_DESCRIPTIONS = {
    "b2": "BCHCode(GF(2),315,128), shortened by 2",
    "ct4": "CodeTables stored generator matrix",
    "b8": "BCHCode(GF(8),189,129,1), shortened by 5",
    "b16": "BCHCode(GF(16),153,128,10)",
    "pl64": "projective-line multiplicity construction",
    "rep": "repetition code",
}


def tex_escape(value: object) -> str:
    return (
        str(value)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("&", r"\&")
        .replace("%", r"\%")
    )


def ensure_dirs() -> None:
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def fetch_codetables_q4_matrix(
    refresh: bool = False, allow_unverified_tls: bool = False
) -> list[list[str]]:
    matrix_path = DATA_DIR / "codetables_q4_177_6_matrix.csv"
    if matrix_path.exists() and not refresh:
        with matrix_path.open("r", newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        validate_q4_matrix(rows)
        return rows

    context = ssl._create_unverified_context() if allow_unverified_tls else None
    try:
        html = urllib.request.urlopen(
            CODETABLES_Q4_URL, timeout=30, context=context
        ).read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError) as exc:
        raise RuntimeError(
            "Could not fetch the GF(4) generator matrix from CodeTables. "
            "Network access and a valid TLS certificate chain are required; "
            "--allow-unverified-tls permits retrieval with the pinned matrix hash."
        ) from exc
    marker = "Construction from a stored generator matrix:"
    if marker not in html:
        raise RuntimeError("CodeTables page did not contain the stored generator matrix marker")
    tail = html.split(marker, 1)[1]
    blocks = re.findall(r"\[\s*(.*?)\s*\]", tail, re.DOTALL)[:6]
    rows = []
    for block in blocks:
        tokens = [token.strip() for token in block.replace("\n", " ").split(",") if token.strip()]
        rows.append(tokens)
    validate_q4_matrix(rows)
    with matrix_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)
    return rows


def validate_q4_matrix(rows: list[list[str]]) -> None:
    allowed = {"0", "1", "w", "w^2"}
    if len(rows) != 6 or any(len(row) != 177 for row in rows):
        raise RuntimeError("Expected a 6 x 177 CodeTables generator matrix")
    bad = sorted({token for row in rows for token in row if token not in allowed})
    if bad:
        raise RuntimeError(f"Unexpected GF(4) entries in CodeTables matrix: {bad}")
    canonical = "\n".join(",".join(row) for row in rows) + "\n"
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    if digest != CODETABLES_Q4_MATRIX_SHA256:
        raise RuntimeError(
            "CodeTables GF(4) matrix differs from the version used in the paper: "
            f"SHA-256 {digest}"
        )


def q4_token_to_magma(token: str) -> str:
    return {"0": "F4!0", "1": "F4!1", "w": "a", "w^2": "a^2"}[token]


def build_magma_script(q4_rows: list[list[str]]) -> str:
    q4_tokens = [q4_token_to_magma(token) for row in q4_rows for token in row]
    q4_entry_lines = []
    for i in range(0, len(q4_tokens), 24):
        q4_entry_lines.append("  " + ", ".join(q4_tokens[i : i + 24]))
    q4_entries = ",\n".join(q4_entry_lines)
    return f"""
SetSeed(1);

procedure Emit(label, q, expected_n, expected_k, target_d, C, method)
  d := MinimumDistance(C);
  ok := (Length(C) eq expected_n) and (Dimension(C) eq expected_k) and (d ge target_d);
  print "CERT|" cat label cat "|" cat Sprint(q) cat "|" cat Sprint(expected_n) cat "|"
        cat Sprint(Length(C)) cat "|" cat Sprint(expected_k) cat "|" cat Sprint(Dimension(C)) cat "|"
        cat Sprint(target_d) cat "|" cat Sprint(d) cat "|" cat Sprint(ok) cat "|" cat method;
end procedure;

F2 := GF(2);
C2parent := BCHCode(F2, 315, 128);
C2 := ShortenCode(C2parent, {{1..2}});
Emit("q2", 2, 313, 12, 128, C2, "b2");

F4<a> := GF(4);
G4 := Matrix(F4, 6, 177, [{q4_entries}]);
C4 := LinearCode(G4);
Emit("q4", 4, 177, 6, 128, C4, "ct4");

F8 := GF(8);
C8parent := BCHCode(F8, 189, 129, 1);
C8 := ShortenCode(C8parent, {{1..5}});
Emit("q8", 8, 184, 4, 129, C8, "b8");

F16 := GF(16);
C16 := BCHCode(F16, 153, 128, 10);
Emit("q16", 16, 153, 3, 128, C16, "b16");

F64 := GF(64);
pts := [];
for a64 in F64 do
  Append(~pts, <F64!1, a64>);
end for;
Append(~pts, <F64!0, F64!1>);
row1 := [];
row2 := [];
for i in [1..#pts] do
  reps := 2;
  if i le 5 then
    reps := 3;
  end if;
  for j in [1..reps] do
    Append(~row1, pts[i][1]);
    Append(~row2, pts[i][2]);
  end for;
end for;
G64 := Matrix(F64, 2, #row1, row1 cat row2);
C64 := LinearCode(G64);
Emit("q64", 64, 135, 2, 132, C64, "pl64");

F4096 := GF(4096);
G4096 := Matrix(F4096, 1, 132, [F4096!1 : i in [1..132]]);
C4096 := LinearCode(G4096);
Emit("q4096", 4096, 132, 1, 132, C4096, "rep");

quit;
"""


def ps_quote(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def run_magma(magma_path: Path, script_path: Path, timeout: int) -> str:
    if os.name == "nt":
        with tempfile.TemporaryDirectory(prefix="qary_ot_magma_") as temp_dir:
            temp_path = Path(temp_dir)
            input_path = temp_path / "input.m"
            output_path = temp_path / "stdout.txt"
            error_path = temp_path / "stderr.txt"
            input_path.write_bytes(script_path.read_bytes())
            command = (
                "$ErrorActionPreference = 'Stop'; "
                f"$p = Start-Process -FilePath {ps_quote(magma_path)} "
                f"-WorkingDirectory {ps_quote(magma_path.parent)} "
                f"-RedirectStandardInput {ps_quote(input_path)} "
                f"-RedirectStandardOutput {ps_quote(output_path)} "
                f"-RedirectStandardError {ps_quote(error_path)} "
                "-WindowStyle Hidden -PassThru -Wait; "
                "exit $p.ExitCode"
            )
            proc = subprocess.run(
                ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", command],
                text=True,
                capture_output=True,
                timeout=timeout,
            )
            output = (
                output_path.read_text(encoding="utf-8", errors="replace")
                if output_path.exists() else ""
            )
            output += (
                error_path.read_text(encoding="utf-8", errors="replace")
                if error_path.exists() else ""
            )
            output += proc.stdout + proc.stderr
    else:
        proc = subprocess.run(
            [str(magma_path)],
            input=script_path.read_text(encoding="utf-8"),
            text=True,
            capture_output=True,
            timeout=timeout,
        )
        output = proc.stdout + proc.stderr
    if proc.returncode != 0:
        raise RuntimeError(f"Magma exited with status {proc.returncode}\n{output}")
    return output


def parse_certificates(transcript: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for line in transcript.splitlines():
        if not line.startswith("CERT|"):
            continue
        parts = line.split("|")
        if len(parts) != 11:
            raise RuntimeError(f"Malformed Magma certificate line: {line}")
        rows.append(
            {
                "label": parts[1],
                "q": parts[2],
                "expected_n": parts[3],
                "magma_n": parts[4],
                "expected_k": parts[5],
                "magma_k": parts[6],
                "target_d": parts[7],
                "magma_min_d": parts[8],
                "ok": parts[9],
                "method": parts[10],
                "method_description": METHOD_DESCRIPTIONS.get(parts[10], parts[10]),
            }
        )
    if not rows:
        raise RuntimeError("Magma transcript did not contain any CERT lines")
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "label",
        "q",
        "expected_n",
        "magma_n",
        "expected_k",
        "magma_k",
        "target_d",
        "magma_min_d",
        "ok",
        "method",
        "method_description",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_latex_table(path: Path, rows: list[dict[str, str]]) -> None:
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\small",
        r"\caption{Magma construction certificates for the selected $\delta=12$ code candidates. The column $d_{\rm Magma}$ is the minimum distance returned by Magma for the constructed code.}",
        r"\label{tab:magma-certificates}",
        r"\begin{tabular}{rrrrrL{.42\linewidth}}",
        r"\toprule",
        r"$q$ & $n$ & $k$ & $d^\star$ & $d_{\rm Magma}$ & source \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            f"{row['q']} & {row['magma_n']} & {row['magma_k']} & {row['target_d']} & "
            f"{row['magma_min_d']} & {tex_escape(row['method_description'])} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Magma certificates for selected q-ary code candidates.")
    parser.add_argument("--magma", type=Path, default=DEFAULT_MAGMA)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--refresh-codetables", action="store_true")
    parser.add_argument(
        "--allow-unverified-tls",
        action="store_true",
        help="Skip TLS certificate validation for CodeTables; the pinned matrix hash is still checked.",
    )
    args = parser.parse_args()

    ensure_dirs()
    magma_executable = str(args.magma)
    resolved_magma = (
        args.magma if args.magma.is_file() else shutil.which(magma_executable)
    )
    if resolved_magma is None:
        raise FileNotFoundError(
            f"Magma executable not found: {magma_executable}; pass --magma PATH or add it to PATH"
        )

    q4_rows = fetch_codetables_q4_matrix(
        refresh=args.refresh_codetables,
        allow_unverified_tls=args.allow_unverified_tls,
    )
    magma_script = build_magma_script(q4_rows)
    input_path = DATA_DIR / "magma_certificates_input.m"
    input_path.write_text(magma_script, encoding="utf-8")
    transcript = run_magma(Path(resolved_magma), input_path, timeout=args.timeout)
    (DATA_DIR / "magma_certificates.log").write_text(transcript, encoding="utf-8")
    rows = parse_certificates(transcript)
    write_csv(DATA_DIR / "magma_certificates.csv", rows)
    write_latex_table(GENERATED_DIR / "magma_certificates_table.tex", rows)

    failed = [row for row in rows if row["ok"] != "true"]
    if failed:
        raise RuntimeError(f"Magma certificates failed: {failed}")

    for row in rows:
        print(
            f"OK q={row['q']} [{row['magma_n']},{row['magma_k']},{row['magma_min_d']}] "
            f"target_d={row['target_d']} {row['method_description']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
