<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/logo/lambda-logo-dark.svg">
    <img src="docs/logo/lambda-logo-light.svg" alt="LAMBDA" width="420">
  </picture>
</p>

# LAMBDA

The executable is `lambdaDrop`.

Eigenvalues of the linear shape oscillations of a free, spherical, viscous drop
(no base flow, no surrounding fluid), computed directly from the linearized
Navier–Stokes equations with a deformable surface.

For each degree *l* the poloidal perturbation is a one-dimensional problem in
the radius, discretized by Chebyshev collocation and solved as a generalized
eigenvalue problem *A x = λ B x*. The solver uses neither the Bessel-function
solution nor Chandrasekhar's characteristic equation, so it gives an
independent check of the roots of that equation:

* relative error below 1e-12 in λ and below 1e-12 in the damping rate Γ for
  1e-3 ≲ Oh ≲ 1 (default `Nr = 48`), see [Validation](#validation);
* 0.01 s per eigenvalue for Oh ≤ 0.02, and 0.1–0.3 s for larger Oh, where the
  root is followed from Oh = 0.01 (measured on two cores).

LAMBDA accompanies

> S. Shiratori et al., "Explicit Corrections for Finite Viscosity and Ambient
> Gas to the Lamb and Rayleigh Formulas for the Oscillating-Drop Method",
> *Int. J. Thermophys.* (submitted).

## Web version (no installation)

`docs/index.html` is a self-contained page that runs the same method in the
browser: choose a material or enter ρ, μ, σ, R (or Oh), the degree *l* and the
order *m*, and the page shows Γ, f and their deviations from Lamb's and
Rayleigh's values, with an animated 3-D view of the mode (rotate with the
mouse, amplitude slider, velocity vectors, a meridional section coloured by
speed, pressure or vorticity). Open the file in a browser, or use the GitHub
Pages site of this repository. The JavaScript solver (`web/src/solver.js`) is
tested against the same exact roots as the C++ code (`cd web && npm test`);
`npm run build` regenerates the page.

## Build

Requirements: a C++20 compiler, CMake ≥ 3.20, [Eigen](https://eigen.tuxfamily.org) 3.4
and [yaml-cpp](https://github.com/jbeder/yaml-cpp) ≥ 0.7. If Eigen or yaml-cpp
is not installed, CMake downloads it (`-DLAMBDA_FETCH_DEPS=OFF` to forbid).
No BLAS/LAPACK is needed.

```sh
cmake -B build
cmake --build build -j
ctest --test-dir build --output-on-failure      # a few seconds
```

* Ubuntu/Debian: `sudo apt install cmake g++ libeigen3-dev libyaml-cpp-dev`
* macOS (Homebrew): `brew install cmake eigen yaml-cpp`
* Windows: Visual Studio 2022 or later, or its Build Tools, with the
  "C++ CMake tools for Windows" component; run the commands above in a
  Developer Command Prompt. The dependencies are downloaded (no git needed),
  or use `vcpkg install eigen3 yaml-cpp`.

## Usage

### YAML input

```sh
./build/lambdaDrop examples/nickel.yaml
```

```yaml
name: Ni_1728K
drop:                        # exactly one of: physical, Oh, We (= Oh^2)
  physical:
    density: 7890            # kg/m^3
    viscosity: 4.7e-3        # Pa s
    surface_tension: 1.739   # N/m
    radius: 1.0e-3           # m
modes: [2, 3, 4]             # degrees l; or [{l: 4, guess: [re, im]}, ...]
resolution:
  Nr: 48                     # Chebyshev functions per variable (default 48)
  check: true                # also solve with Nr + 16, report the changes of Γ and ω
solver:
  method: mode               # mode (default) | spectrum (all eigenvalues)
  tolerance: 1.0e-10         # relative backward error to accept
output:
  file: nickel_out.yaml      # .yaml or .csv
  eigenfunction:             # optional: spatial form of the modes
    m: [0, 2]                # azimuthal orders 0 <= m <= l
    format: [ascii, vtk]     # .dat tables and/or ParaView .vts
```

The output gives λ = −Γ + iω in viscous units (time R²/ν) and capillary units
(time (ρR³/σ)^½), and, for physical input, Γ in s⁻¹ and the frequency
f = ω/2π in Hz. Unknown keys are rejected.

### Eigenfunctions

With `output.eigenfunction` the spatial form of each mode is written for each
requested *m* (see `examples/eigenfunction.yaml`):

| file | content |
|---|---|
| `<prefix>_l<l>_radial.dat` | u(r), v(r), p(r) (u_r = uY, u_tan = v r∇Y, p = pY) |
| `<prefix>_l<l>_m<m>_meridional.dat` | velocity and pressure in the plane y = 0 |
| `<prefix>_l<l>_m<m>_surface.dat` | surface displacement dY and surface velocity |
| `<prefix>_l<l>_m<m>.vts` | velocity and pressure in the drop (ParaView) |
| `<prefix>_l<l>_m<m>_surface.vts` | sphere with the displacement d Y e_r (ParaView) |

The amplitudes are complex and normalised so that the surface is
r = 1 + Re[d Y e^{λt}] with d = 1 and max|Y| = 1, Y = P_l^m(cos θ) e^{imφ}
(Condon–Shortley phase); every quantity q is given as `_re` and `_im`, and the
field at time t is e^{−Γt}(Re q cos ωt − Im q sin ωt). Units: length R,
velocity ν/R, pressure ρν²/R². In ParaView, apply *Warp By Vector* with
`displacement_vector_re` and a small scale factor to see the deformed surface;
use the Calculator to combine `_re` and `_im` for other phases. Output
resolution: `points: {r: 33, theta: 49, phi: 65}` (default).

### Command line

```sh
./build/lambdaDrop We l [m] [Nr] [guessRe guessIm]   # one eigenvalue
./build/lambdaDrop --all We l [Nr]                    # all eigenvalues of degree l
```

`We = Oh² = μ²/(ρσR)`; λ is printed in viscous units. The azimuthal order
`m` may be given but does not change the eigenvalue of a sphere.

### Which root is returned

Without a guess, the shape mode of degree *l* is found from Lamb's damping
rate and Rayleigh's frequency for Oh ≤ 0.02, and for larger Oh by following
the root from Oh = 0.01. Past the critical Oh the oscillation becomes
aperiodic and the pair turns into two real roots; the slower one (least
damped) is returned. With a guess, the eigenvalue nearest to it is returned.

### Resolution

The boundary layer at the surface has thickness ~ (ν/ω)^½, so small Oh
and large *l* need more radial functions. With `Nr = 32` the error in Γ is
1.6e-9 for Oh = 1.27e-3, *l* = 4; with `Nr = 48` it is 3e-13. Use
`resolution.check: true` to see the relative changes of Γ and of ω with
`Nr + 16` (reported separately, since Γ ≪ ω for small Oh).

## Method

See [docs/formulation.md](docs/formulation.md). In short: velocity
u = u(r) Y e_r + v(r) r∇Y, pressure p(r) Y and surface r = 1 + d Y
(Y = Y_l^m); continuity and the radial and tangential momentum equations are
collocated at Gauss–Lobatto points in r ∈ (0, 1]; at r = 1 the momentum
equations are replaced by the normal-stress and tangential-stress conditions,
and the kinematic condition gives the equation for d. Regularity at r = 0 is
built into the basis (Chebyshev polynomials of the parity of r^(l−1) for u, v
and of r^l for p). All finite eigenvalues are computed by the QZ algorithm and
the selected one is refined by shifted inverse iteration, the shift being
updated at every step by a Rayleigh quotient. The results are deterministic. For a sphere the
eigenvalues do not depend on the azimuthal order m.

## Validation

`tests/reference/chandrasekhar_roots.csv` holds 25 roots of Chandrasekhar's
dispersion relation computed with 40 digits (Python, mpmath; *l* = 2–4,
Oh = 1.26e-3 to 1). It is generated by `tools/make_reference.py`, which uses
the root finder in `tools/chandrasekhar_roots.py`; CI checks that it is
reproduced. The test `chandrasekhar` requires relative errors below 1e-12 in
λ, Γ and ω at `Nr = 48` (the largest is 4.3e-13, in Γ); the test
`convergence` checks the spectral convergence in `Nr`; the test `eigenfunction` checks the Legendre
functions, the kinematic and stress-free conditions, div u = 0 and the
potential-flow limit of the eigenfunctions.

```sh
python3 tools/chandrasekhar_roots.py --oh 0.0316 --l 2 3 4    # needs mpmath
```

## Citation

Software: S. Shiratori, *LAMBDA* (version 1.0.0), Zenodo, DOI to be added at
the release; see [CITATION.cff](CITATION.cff). Please also cite the article
given at the top of this file.

## Licence

MIT, see [LICENSE](LICENSE).
