// Spectral convergence in Nr: the error must fall from Nr = 24 to 32 to 48 for the
// thinnest boundary layer of the validation set (Ni, Oh = 1.27e-3, l = 4), and stay
// below 1e-12 at Nr = 48 and 64. At Nr = 96 rounding errors grow (about 2e-12);
// the bound there is 1e-11.
#include <cmath>
#include <cstdio>

#include "solver/Solver.h"

using namespace lambda;

int main() {
    const double Oh = std::sqrt(1.60997e-6);
    const cplx ex(-26.21961135891691720091, 6686.610778964921973027);  // reference CSV
    int fail = 0;
    double prev = 1.0;
    for (int Nr : {24, 32, 48, 64, 96}) {
        const auto r = solve_mode(4, Oh, Nr, std::nullopt, 1e-10);
        const double eG = std::abs(r.lambda.real() / ex.real() - 1.0);
        std::printf("Nr %3d  dGamma/Gamma %.2e\n", Nr, eG);
        if (Nr <= 48 && !(eG < prev)) ++fail;  // decreasing
        if ((Nr == 48 || Nr == 64) && !(eG < 1e-12)) ++fail;  // converged
        if (Nr == 96 && !(eG < 1e-11)) ++fail;                // rounding
        prev = eG;
    }
    std::printf(fail ? "FAIL\n" : "ok\n");
    return fail ? 1 : 0;
}
