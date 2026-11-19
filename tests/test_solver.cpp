// Compares LAMBDA with exact roots of Chandrasekhar's dispersion relation
// (tests/reference/chandrasekhar_roots.csv, made by tools/make_reference.py).
// Criterion at Nr = 48: relative errors of lambda, Gamma and omega below 1e-12.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "solver/Solver.h"

using namespace lambda;

int main() {
    std::ifstream in(LAMBDA_REFERENCE_CSV);
    if (!in) {
        std::printf("cannot open %s\n", LAMBDA_REFERENCE_CSV);
        return 1;
    }
    constexpr int Nr = 48;
    constexpr double tol_lambda = 1e-12;  // |dlambda|/|lambda|
    constexpr double tol_gamma = 1e-12;   // |dGamma|/Gamma
    int n = 0, nfail = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'W') continue;
        std::stringstream ss(line);
        std::string w, l, re, im;
        std::getline(ss, w, ',');
        std::getline(ss, l, ',');
        std::getline(ss, re, ',');
        std::getline(ss, im, ',');
        const double We = std::stod(w);
        const cplx ex(std::stod(re), std::stod(im));
        const ModeResult r = solve_mode(std::stoi(l), std::sqrt(We), Nr, std::nullopt, 1e-10);
        const double eL = std::abs(r.lambda - ex) / std::abs(ex);
        const double eG = std::abs(r.lambda.real() / ex.real() - 1.0);
        const double eW = ex.imag() != 0 ? std::abs(r.lambda.imag() / ex.imag() - 1.0)
                                         : std::abs(r.lambda.imag()) / std::abs(ex);
        const bool ok = r.converged && eL < tol_lambda && eG < tol_gamma && eW < tol_gamma;
        std::printf(
            "We %-10s l %s  lambda %+.15e %+.15e  |dl|/|l| %.1e  dG/G %.1e  dw/w %.1e  bwd %.1e  %s\n",
            w.c_str(), l.c_str(), r.lambda.real(), r.lambda.imag(), eL, eG, eW, r.backward_error,
            ok ? "ok" : "FAIL");
        ++n;
        if (!ok) ++nfail;
    }
    std::printf("%d of %d passed\n", n - nfail, n);
    return (n > 0 && nfail == 0) ? 0 : 1;
}
