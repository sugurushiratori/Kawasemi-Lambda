// Linear shape oscillations of a free viscous drop: for each degree l the
// poloidal perturbation reduces to a one-dimensional problem in r, discretised
// by Chebyshev collocation (see docs/formulation.md). The eigenvalues do not
// depend on the azimuthal order m (spherical symmetry).
#pragma once

#include <complex>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include "io/Config.h"

namespace lambda {

using cplx = std::complex<double>;

struct ModeResult {
    int l = 0;
    cplx lambda;                  // viscous units (time R^2/nu); -Gamma + i omega
    double backward_error = 0.0;  // ||A x - lambda B x|| / (||A x|| + |lambda| ||B x||)
    int iterations = 0;           // refinement steps
    bool converged = false;
    // relative changes with Nr -> Nr + dNr, if checked (< 0: not computed;
    // the omega change is not computed for an aperiodic root, omega = 0)
    double resolution_change_gamma = -1.0;  // |Gamma(Nr+dNr) - Gamma(Nr)| / Gamma(Nr)
    double resolution_change_omega = -1.0;  // |omega(Nr+dNr) - omega(Nr)| / omega(Nr)

    // Eigenvector: Chebyshev coefficients [u (Nr), v (Nr), p (Nr), d], normalised to d = 1.
    int Nr = 0;
    Eigen::VectorXcd coeffs;
};

// Collocation matrices A x = lambda B x (rows scaled to unit max norm).
struct Pencil {
    Eigen::MatrixXd A, B;
};
Pencil poloidal_pencil(int l, double Oh, int Nr);

// Lamb's damping with Rayleigh's frequency, in viscous units (default guess).
cplx lamb_guess(int l, double Oh);

// All finite eigenvalues of the l-th poloidal problem, sorted by decreasing real part.
std::vector<cplx> spectrum(int l, double Oh, int Nr);

// The shape-oscillation eigenvalue of degree l.
//   guess given: the eigenvalue nearest to it;
//   otherwise  : nearest to Lamb's value for Oh <= 0.02, or followed in Oh
//                from Oh = 0.01 for larger Oh. Beyond the critical Oh the pair
//                becomes two real roots; the slower (least damped) one is returned.
ModeResult solve_mode(int l, double Oh, int Nr, const std::optional<cplx>& guess,
                      double tolerance);

std::vector<ModeResult> solve(const Config& cfg);

}  // namespace lambda
