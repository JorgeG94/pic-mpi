#!/usr/bin/env python3
"""Require the serial backend to present the same procedures as mpi_f08.

    python3 tools/check_serial_in_step.py

WHY THIS IS CHECKED RATHER THAN REMEMBERED
------------------------------------------
`src/serial/pic_mpi_serial.F90` is generated from `src/mpi_f08/pic_mpi.F90`
by `tools/gen_serial.py`: same 140 signatures, different bodies. Generated
code that lives in the repository drifts from its generator the moment someone
edits one and not the other.

The failure mode is nasty and remote. Add a procedure to the MPI backend, do
not regenerate, and nothing here breaks -- every MPI build is fine. It breaks
in a CONSUMER's serial build, as a generic that will not resolve, with no clue
pointing back at this repository.

So the two are compared here: same procedure names, same dummy argument lists.
Bodies are expected to differ and are not compared.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

PROC_RE = re.compile(r'^\s*(pure\s+|elemental\s+)?(subroutine|function)\s+(\w+)\s*\(([^)]*)',
                     re.I)


def signatures(path: Path) -> dict[str, list[str]]:
    """{name: [dummy, ...]} for every procedure after the module's contains."""
    lines = path.read_text().splitlines()
    depth, ci = 0, None
    for i, l in enumerate(lines):
        s = l.strip().lower()
        if re.match(r'type\s*(,.*)?::\s*\w+', s) and 'end type' not in s:
            depth += 1
        elif s.startswith('end type'):
            depth -= 1
        elif s == 'contains' and depth == 0:
            ci = i
            break
    if ci is None:
        sys.exit(f"{path}: no module-level `contains` found")

    out: dict[str, list[str]] = {}
    i = ci + 1
    while i < len(lines):
        m = PROC_RE.match(lines[i])
        if m:
            # join continuations so the argument list is complete
            txt, j = lines[i], i
            while txt.rstrip().endswith('&'):
                j += 1
                txt = txt.rstrip()[:-1] + lines[j]
            am = re.search(r'\((.*?)\)', txt)
            args = [a.strip().lower() for a in am.group(1).split(',')] if am else []
            out[m.group(3).lower()] = args
            i = j
        i += 1
    return out


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    mpi = signatures(root / 'src' / 'mpi_f08' / 'pic_mpi.F90')
    ser = signatures(root / 'src' / 'serial' / 'pic_mpi_serial.F90')

    problems: list[str] = []
    for name in sorted(set(mpi) - set(ser)):
        problems.append(f"missing from the serial backend: {name}")
    for name in sorted(set(ser) - set(mpi)):
        problems.append(f"only in the serial backend: {name}")
    for name in sorted(set(mpi) & set(ser)):
        if mpi[name] != ser[name]:
            problems.append(
                f"{name}: arguments differ\n"
                f"      mpi_f08 ({len(mpi[name])}): {', '.join(mpi[name])}\n"
                f"      serial  ({len(ser[name])}): {', '.join(ser[name])}")

    print(f"  mpi_f08: {len(mpi)} procedures")
    print(f"  serial : {len(ser)} procedures")
    for p in problems:
        print(f"  {p}")
    if problems:
        print(f"\n  {len(problems)} difference(s). Regenerate with:")
        print("    python3 tools/gen_serial.py src/mpi_f08/pic_mpi.F90 \\")
        print("            src/serial/pic_mpi_serial.F90")
        return 1
    print("\n  the two backends present the same interface")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
