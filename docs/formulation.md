# Formulation

## Problem

A drop of radius R, density ρ, viscosity μ (ν = μ/ρ) and surface tension σ,
in vacuum and at rest. Small perturbations ∝ exp(λt) of the velocity, the
pressure and the surface position r = R(1 + d Y_l^m(θ, φ)).

Units: length R, time R²/ν, pressure ρν²/R². The only parameter is the
Ohnesorge number Oh = μ/(ρσR)^½; the command line uses We = Oh².
λ = −Γ + iω. Capillary units (time (ρR³/σ)^½) follow by multiplying λ by Oh.

## Reduction to one dimension

For a sphere the linear problem separates in spherical harmonics. The
poloidal part of degree l is

    u = u(r) Y e_r + v(r) r ∇Y,      p = p(r) Y,      L = l(l+1).

(The toroidal part decouples, does not move the surface and is not computed.)

Equations, for 0 < r < 1:

    continuity        u' + 2u/r − L v/r                    = 0
    radial momentum   −L u/r² + L (v' + v/r)/r − p'         = λ u
    tangential mom.   −u'/r + v'' + 2v'/r − p/r             = λ v

The viscous terms use ∇²u = −∇×∇×u for a solenoidal field.

Conditions at r = 1:

    normal stress     p − 2u' − (l−1)(l+2) d / Oh²         = 0
    tangential stress v' − v + u                           = 0
    kinematic         u                                    = λ d

## Discretization

Regularity at r = 0 requires u, v ∝ r^(l−1) and p ∝ r^l. Each function is
expanded in Chebyshev polynomials of one parity on [−1, 1],

    u = Σ a_n T_{2n+q}(r),   n = 0 … Nr−1,   q = (l+1) mod 2   (u, v)
                                             q = l mod 2       (p)

so that the expansion has the parity of r^(l−1) or r^l.

Collocation points: r_j = cos(π/2 · (1 − j/(Nr−1))), j = 1 … Nr−1 (the
Gauss–Lobatto points of [−1, 1] with r > 0, r_{Nr−1} = 1), and
r_0 = cos(π/2 · (1 − 1/(2(Nr−1)))) in place of r = 0.

Rows: continuity at all Nr points; the two momentum equations at the Nr−1
interior points; the two stress conditions at r = 1; the kinematic condition.
Unknowns: 3Nr coefficients and d, so A and B are square of size 3Nr+1.
B is singular (continuity and the boundary conditions contain no λ), which
gives infinite eigenvalues; each row of A and B is divided by its largest
entry.

## Eigenvalues

All finite eigenvalues are computed with the real QZ algorithm (Eigen's
`GeneralizedEigenSolver`). The selected one is refined in complex arithmetic
by shifted inverse iteration, (A − σB) y = B x, x ← y/‖y‖, with the shift
updated at every step by the one-sided (least-squares) Rayleigh quotient
σ = (Bx)ᴴAx / ‖Bx‖². No random numbers are used, so the results are
deterministic. The result is accepted when the relative backward error

    η = ‖A x − λ B x‖ / (‖A x‖ + |λ| ‖B x‖)

is below `tolerance` (default 1e-10; typically 1e-16 is reached).

Root selection: see README, "Which root is returned".

## Eigenfunctions

The eigenvector (Chebyshev coefficients of u, v, p and d) is divided by d,
so the surface displacement amplitude is 1. The angular factor is
Y = P_l^m(cos θ) e^{imφ} / max_θ |P_l^m| (Condon–Shortley phase), so that
max |Y| = 1; then

    u_r = u(r) Y,   u_θ = v(r) ∂Y/∂θ,   u_φ = i m v(r) Y / sin θ,   p = p(r) Y.

dP_l^m/dθ and P_l^m/sin θ are computed from the derivatives of P_l in cos θ,
so they are finite at the poles.

## Reference: the exact dispersion relation

`tools/chandrasekhar_roots.py` solves, with q² = λ and
Q_l(q) = q i_{l+1}(q)/i_l(q) (modified spherical Bessel functions),

    ω_R² (2Q − q²) = q² [ q⁴ + 2(l−1)(2l+1) q² − 2Q q² − 4 l(l−1)(l+2) Q ],
    ω_R² = l(l−1)(l+2)/Oh²,

which is the free-drop limit of Chandrasekhar's characteristic equation
(Chandrasekhar 1959; Reid 1960). Q_l is even in q, so there is no branch cut.

* S. Chandrasekhar, The oscillations of a viscous liquid globe,
  Proc. London Math. Soc. (3) 9 (1959) 141–149.
* W. H. Reid, The oscillations of a viscous liquid drop,
  Q. Appl. Math. 18 (1960) 86–89.
