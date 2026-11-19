#!/usr/bin/env python3
"""Write tests/reference/chandrasekhar_roots.csv: exact roots for the solver tests.

Cases: the six drops of the validation table (Pt, Ni, Al2O3, Oh = 0.0316, slag,
Oh = 0.316; parameter We = Oh^2) for l = 2, 3, 4, plus large Oh where the roots
are followed from small Oh (including aperiodic roots past the critical Oh).
"""
import os, sys
import mpmath as mp
sys.path.insert(0, os.path.dirname(__file__))
from chandrasekhar_roots import root, tracked_root

mp.mp.dps = 40
WE = ["1.60997e-6", "1.5884e-6", "4.61538e-4", "1e-3", "4.0386e-2", "1e-1"]
EXTRA = [("0.25", [2, 3, 4]), ("0.49", [2, 3, 4]), ("1", [2])]

out = os.path.join(os.path.dirname(__file__), "..", "tests", "reference", "chandrasekhar_roots.csv")
with open(out, "w") as f:
    f.write("# Exact roots of Chandrasekhar's dispersion relation (tools/make_reference.py, mpmath 40 digits)\n")
    f.write("# We = Oh^2; lambda = -Gamma + i omega in viscous units (time R^2/nu)\n")
    f.write("We,l,re,im\n")
    rows = [(w, l) for w in WE for l in (2, 3, 4)] + [(w, l) for w, ls in EXTRA for l in ls]
    for w, l in rows:
        oh = mp.sqrt(mp.mpf(w))
        lam = tracked_root(l, oh) if oh > 0.02 else root(l, oh)
        if abs(lam.imag) < mp.mpf(10) ** -20 * abs(lam):
            lam = mp.mpc(lam.real, 0)
        f.write(f"{w},{l},{mp.nstr(lam.real, 22)},{mp.nstr(lam.imag, 22)}\n")
print("wrote", os.path.normpath(out))
