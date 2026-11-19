#include "solver/Solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "solver/RadialBasis.h"

namespace lambda {

// ---------------------------------------------------------------------------
//  Poloidal problem of degree l. Units: length R, time R^2/nu, pressure
//  rho nu^2/R^2. Perturbation ~ exp(lambda t) with
//      u_r = u(r) Y_l^m,   u_tan = v(r) r grad Y_l^m,   p = p(r) Y_l^m,
//      surface r = 1 + d Y_l^m,           L = l(l+1),   Oh^2 = mu^2/(rho sigma R).
//  Interior points r_j (j < Nr-1):
//      continuity     u' + 2u/r - L v/r                          = 0
//      r momentum     -L u/r^2 + L (v' + v/r)/r - p'              = lambda u
//      tang. momentum -u'/r + v'' + 2v'/r - p/r                   = lambda v
//  Surface r = 1 (replacing the two momentum rows):
//      normal stress  p - 2u' - (l-1)(l+2) d / Oh^2               = 0
//      tang. stress   v' - v + u                                  = 0
//      kinematic      u(1)                                        = lambda d
//  Continuity is imposed at every point, including r = 1.
// ---------------------------------------------------------------------------
Pencil poloidal_pencil(int l, double Oh, int Nr) {
    if (l < 2) throw std::invalid_argument("l must be >= 2");
    const Eigen::VectorXd r = RadialBasis::points(Nr);
    const RadialBasis V(Nr, (l + 1) % 2);  // u, v ~ r^(l-1)
    const RadialBasis S(Nr, l % 2);        // p    ~ r^l
    const double L = l * (l + 1.0);
    const double We = Oh * Oh;

    const int n = 3 * Nr + 1;
    const int iu = 0, iv = Nr, ip = 2 * Nr, id = 3 * Nr;
    Pencil P{Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n)};
    auto& A = P.A;
    auto& B = P.B;

    for (int j = 0; j < Nr; ++j) {
        const double x = r[j];
        const bool surface = (j == Nr - 1);
        const int rc = j, rr = Nr + j, rt = 2 * Nr + j;
        for (int k = 0; k < Nr; ++k) {
            const double u0 = V.T(k, j), u1 = V.T1(k, j), u2 = V.T2(k, j);
            const double p0 = S.T(k, j), p1 = S.T1(k, j);
            // continuity
            A(rc, iu + k) = u1 + 2.0 * u0 / x;
            A(rc, iv + k) = -L * u0 / x;
            if (!surface) {
                // r momentum
                A(rr, iu + k) = -L / (x * x) * u0;
                B(rr, iu + k) = u0;
                A(rr, iv + k) = L / x * (u1 + u0 / x);
                A(rr, ip + k) = -p1;
                // tangential momentum
                A(rt, iu + k) = -u1 / x;
                A(rt, iv + k) = u2 + 2.0 * u1 / x;
                B(rt, iv + k) = u0;
                A(rt, ip + k) = -p0 / x;
            } else {
                A(rr, ip + k) = p0;
                A(rr, iu + k) = -2.0 * u1;
                A(rt, iv + k) = u1 - u0;
                A(rt, iu + k) = u0;
                A(id, iu + k) = u0;  // kinematic
            }
        }
        if (surface) A(rr, id) = -(l - 1.0) * (l + 2.0) / We;
    }
    B(id, id) = 1.0;

    // scale each row by its largest entry in A or B
    for (int i = 0; i < n; ++i) {
        double s = std::max(A.row(i).cwiseAbs().maxCoeff(), B.row(i).cwiseAbs().maxCoeff());
        if (s > 0) {
            A.row(i) /= s;
            B.row(i) /= s;
        }
    }
    return P;
}

cplx lamb_guess(int l, double Oh) {
    return {-(l - 1.0) * (2.0 * l + 1.0), std::sqrt(l * (l - 1.0) * (l + 2.0)) / Oh};
}

namespace {

std::vector<cplx> finite_eigenvalues(const Pencil& P) {
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(P.A, P.B, false);
    if (ges.info() != Eigen::Success) throw std::runtime_error("QZ did not converge");
    const auto alpha = ges.alphas();
    const auto beta = ges.betas();
    const double amax = alpha.cwiseAbs().maxCoeff();
    std::vector<cplx> ev;
    for (Eigen::Index i = 0; i < alpha.size(); ++i) {
        if (std::abs(beta[i]) <= 1e-13 * amax) continue;  // infinite eigenvalue
        cplx z = alpha[i] / beta[i];
        if (std::isfinite(z.real()) && std::isfinite(z.imag())) ev.push_back(z);
    }
    std::sort(ev.begin(), ev.end(), [](cplx a, cplx b) { return a.real() > b.real(); });
    return ev;
}

bool is_real(cplx z) { return std::abs(z.imag()) <= 1e-10 * std::abs(z); }

// Once the oscillatory pair has split into two real roots (aperiodic decay),
// the shape mode is taken as the slower of the two: the least-damped real root.
cplx least_damped_real(const std::vector<cplx>& ev) {
    for (cplx z : ev)  // sorted by decreasing real part
        if (is_real(z)) return {z.real(), 0.0};
    throw std::runtime_error("no real eigenvalue");
}

cplx nearest(const std::vector<cplx>& ev, cplx target) {
    if (ev.empty()) throw std::runtime_error("no finite eigenvalue");
    return *std::min_element(ev.begin(), ev.end(), [&](cplx a, cplx b) {
        return std::abs(a - target) < std::abs(b - target);
    });
}

// Shifted inverse iteration started at sigma, the shift updated at every step by a
// one-sided (least-squares) Rayleigh quotient (Bx)^H A x / ||Bx||^2. Refines the
// eigenvalue from QZ and gives the relative backward error of the final pair.
ModeResult refine(const Pencil& P, cplx sigma, double tolerance) {
    using CM = Eigen::MatrixXcd;
    using CV = Eigen::VectorXcd;
    const CM A = P.A.cast<cplx>(), B = P.B.cast<cplx>();
    const Eigen::Index n = A.rows();
    CV x = CV::Ones(n);
    ModeResult res;
    cplx lam = sigma;
    CV best = x;
    double eta = std::numeric_limits<double>::infinity();
    for (int it = 1; it <= 8; ++it) {
        // a tiny offset keeps A - sigma B invertible when sigma is already exact
        const cplx shift = lam * cplx(1.0 + 1e-14, 0.0);
        Eigen::PartialPivLU<CM> lu(A - shift * B);
        CV y = lu.solve(B * x);
        x = y / y.norm();
        const CV Ax = A * x, Bx = B * x;
        const cplx lam_new = Bx.dot(Ax) / Bx.squaredNorm();  // least-squares eigenvalue
        const double eta_new =
            (Ax - lam_new * Bx).norm() / (Ax.norm() + std::abs(lam_new) * Bx.norm());
        res.iterations = it;
        const bool stalled = eta_new >= 0.5 * eta && it >= 2;
        if (eta_new < eta) {
            lam = lam_new;
            eta = eta_new;
            best = x;
        }
        if (eta < 1e-15 || stalled) break;
    }
    res.lambda = lam;
    res.backward_error = eta;
    const cplx dsurf = best[n - 1];  // surface displacement
    res.coeffs = (std::abs(dsurf) > 0) ? CV(best / dsurf) : best;
    res.converged = eta < tolerance;
    return res;
}

}  // namespace

std::vector<cplx> spectrum(int l, double Oh, int Nr) {
    return finite_eigenvalues(poloidal_pencil(l, Oh, Nr));
}

ModeResult solve_mode(int l, double Oh, int Nr, const std::optional<cplx>& guess,
                      double tolerance) {
    constexpr double kTrackFrom = 0.01;  // Lamb's value is reliable below this Oh
    cplx target;
    if (guess) {
        target = *guess;
    } else if (Oh <= 2.0 * kTrackFrom) {
        target = lamb_guess(l, Oh);
    } else {
        // follow the root in log(Oh) from kTrackFrom
        const int steps = static_cast<int>(std::ceil(20.0 * std::log10(Oh / kTrackFrom)));
        double o = kTrackFrom;
        target = nearest(spectrum(l, o, Nr), lamb_guess(l, o));
        for (int s = 1; s <= steps; ++s) {
            const double o_new = kTrackFrom * std::pow(Oh / kTrackFrom, double(s) / steps);
            // lambda/Oh-scaling of the frequency keeps the prediction close
            const cplx pred(target.real(), target.imag() * o / o_new);
            const auto ev = spectrum(l, o_new, Nr);
            target = is_real(target) ? least_damped_real(ev) : nearest(ev, pred);
            o = o_new;
        }
    }
    const Pencil P = poloidal_pencil(l, Oh, Nr);
    const auto ev = finite_eigenvalues(P);
    const cplx qz = (!guess && is_real(target)) ? least_damped_real(ev) : nearest(ev, target);
    ModeResult res = refine(P, qz, tolerance);
    res.l = l;
    res.Nr = Nr;
    if (res.lambda.imag() < 0) {  // report the member of the pair with omega > 0
        res.lambda = std::conj(res.lambda);
        res.coeffs = res.coeffs.conjugate();
    }
    if (is_real(res.lambda)) res.lambda.imag(0.0);
    return res;
}

std::vector<ModeResult> solve(const Config& cfg) {
    std::vector<ModeResult> out;
    for (const auto& md : cfg.modes) {
        ModeResult r = solve_mode(md.l, cfg.Oh, cfg.Nr, md.guess, cfg.tolerance);
        if (cfg.check_resolution) {
            const ModeResult r2 =
                solve_mode(md.l, cfg.Oh, cfg.Nr + cfg.check_dNr, r.lambda, cfg.tolerance);
            r.resolution_change_gamma =
                std::abs(r2.lambda.real() - r.lambda.real()) / std::abs(r.lambda.real());
            if (r.lambda.imag() > 0)
                r.resolution_change_omega =
                    std::abs(r2.lambda.imag() - r.lambda.imag()) / r.lambda.imag();
        }
        out.push_back(r);
    }
    return out;
}

}  // namespace lambda
