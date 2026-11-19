#!/usr/bin/env python3
"""Exact eigenvalues of a viscous drop oscillating in vacuum (Chandrasekhar 1959, Reid 1960).

Independent reference for the LAMBDA solver. Pure Python + mpmath.

Units: length R, time R^2/nu (viscous units).  lambda = -Gamma + i*omega.
The Ohnesorge number is Oh = mu / sqrt(rho*sigma*R); the solver input "We" is Oh^2.

Dispersion relation, with q^2 = lambda, Q_l(q) = q i_{l+1}(q)/i_l(q),
w2 = l(l-1)(l+2)/Oh^2:

    w2 (2Q - q^2) = q^2 [ q^4 + 2(l-1)(2l+1) q^2 - 2Q q^2 - 4 l(l-1)(l+2) Q ]

Q_l is even in q, so the relation has no branch cut in lambda.

Usage:
    chandrasekhar_roots.py --oh 0.0316 --l 2 3 4
    chandrasekhar_roots.py --we 1e-3 --l 2
    chandrasekhar_roots.py --oh 0.2 --l 4 --no-track   # do not follow the root from small Oh
"""
import argparse
import mpmath as mp


def Q(l, q):
    # q i_{l+1}(q)/i_l(q) = q I_{l+3/2}(q)/I_{l+1/2}(q)
    return q * mp.besseli(l + mp.mpf(3) / 2, q) / mp.besseli(l + mp.mpf(1) / 2, q)


def residual(lam, l, oh):
    q = mp.sqrt(lam)
    Ql = Q(l, q)
    w2 = l * (l - 1) * (l + 2) / oh**2
    q2 = lam
    rhs = q2 * (q2**2 + 2 * (l - 1) * (2 * l + 1) * q2 - 2 * Ql * q2
                - 4 * l * (l - 1) * (l + 2) * Ql)
    return w2 * (2 * Ql - q2) - rhs


def lamb(l, oh):
    """Lamb's weak-damping eigenvalue in viscous units."""
    g = mp.mpf((l - 1) * (2 * l + 1))
    w2 = mp.mpf(l * (l - 1) * (l + 2)) / oh**2
    return mp.mpc(-g, mp.sqrt(max(w2 - g**2, mp.mpf(0))))


def root(l, oh, guess=None):
    if guess is None:
        guess = lamb(l, oh)
    return mp.findroot(lambda z: residual(z, l, oh), mp.mpc(guess), tol=mp.mpf(10)**(-mp.mp.dps + 5))


def tracked_root(l, oh, oh_start=1e-3, steps_per_decade=40):
    """Follow the l-th branch from small Oh (where Lamb is a good guess) up to oh."""
    oh = mp.mpf(oh)
    o = min(mp.mpf(oh_start), oh)
    lam = root(l, o)
    if o == oh:
        return lam
    n = int(mp.ceil(steps_per_decade * mp.log10(oh / o)))
    for k in range(1, n + 1):
        o_new = o * (oh / o) ** (mp.mpf(1) / (n - k + 1))
        # rescale previous root: omega ~ 1/Oh, Gamma ~ const for small Oh
        lam = root(l, o_new, guess=lam)
        o = o_new
    return lam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--oh", type=str, help="Ohnesorge number")
    g.add_argument("--we", type=str, help="solver parameter We = Oh^2")
    ap.add_argument("--l", type=int, nargs="+", default=[2])
    ap.add_argument("--dps", type=int, default=40, help="mpmath decimal digits")
    ap.add_argument("--no-track", action="store_true",
                    help="start from Lamb's value directly; by default the root is followed "
                         "from Oh=1e-3 when Oh > 0.02, because Lamb's guess can land on another branch")
    ap.add_argument("--guess", type=str, nargs=2, metavar=("RE", "IM"))
    a = ap.parse_args()

    mp.mp.dps = a.dps
    oh = mp.mpf(a.oh) if a.oh else mp.sqrt(mp.mpf(a.we))
    print(f"# Oh = {mp.nstr(oh, 15)}  We = {mp.nstr(oh**2, 15)}  (viscous units R^2/nu)")
    print("# l  Re(lambda)                 Im(lambda)                 |residual|")
    for l in a.l:
        if a.guess:
            lam = root(l, oh, mp.mpc(mp.mpf(a.guess[0]), mp.mpf(a.guess[1])))
        elif oh > 0.02 and not a.no_track:
            lam = tracked_root(l, oh)
        else:
            lam = root(l, oh)
        if abs(lam) < 1e-6:
            raise SystemExit(f"l={l}: converged to the trivial root lambda=0; give --guess")
        r = abs(residual(lam, l, oh)) / abs(lam) ** 3
        print(f"{l:3d}  {mp.nstr(lam.real, 20, min_fixed=-30, max_fixed=30):>25s}  "
              f"{mp.nstr(lam.imag, 20, min_fixed=-30, max_fixed=30):>25s}  {mp.nstr(r, 3)}")


if __name__ == "__main__":
    main()
