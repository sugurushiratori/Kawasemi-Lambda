#!/usr/bin/env python3
"""Validation table: LAMBDA against the exact roots (tests/reference/chandrasekhar_roots.csv).

    python3 tools/validation_table.py [path/to/lambdaDrop] [Nr]

Prints We, l, LAMBDA's lambda, the exact root and the relative differences in
Gamma and omega. Needs only the Python standard library.
"""
import csv, os, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
exe = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "build", "lambdaDrop")
Nr = sys.argv[2] if len(sys.argv) > 2 else "48"
ref = os.path.join(here, "..", "tests", "reference", "chandrasekhar_roots.csv")

rows = [r for r in csv.reader(l for l in open(ref) if not l.startswith("#"))][1:]
print(f"# LAMBDA (Nr = {Nr}) vs exact roots; lambda = -Gamma + i omega, viscous units")
print(f"# {'We':>10s} {'l':>2s} {'Re(lambda)':>24s} {'Im(lambda)':>24s} "
      f"{'Re exact':>24s} {'Im exact':>24s} {'dGamma/Gamma':>12s} {'domega/omega':>12s}")
for we, l, re, im in rows:
    out = subprocess.run([exe, we, l, "0", Nr], capture_output=True, text=True, check=True).stdout
    lr, li = map(float, out.strip().splitlines()[-1].split()[1:3])
    ere, eim = float(re), float(im)
    dg = lr / ere - 1.0
    dw = f"{li / eim - 1.0:12.1e}" if eim != 0 else f"{'(real)':>12s}"
    print(f"  {we:>10s} {l:>2s} {lr:24.16e} {li:24.16e} {re:>24.22s} {im:>24.22s} {dg:12.1e} {dw}")
