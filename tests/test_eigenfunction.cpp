// Checks of the eigenfunctions: Legendre functions, kinematic condition,
// potential-flow limit, solenoidality and the stress-free surface.
#include <cmath>
#include <cstdio>
#include <numbers>

#include "post/Eigenfunction.h"
#include "solver/Solver.h"

using namespace lambda;

static int failures = 0;
static void check(bool ok, const char* what, double val) {
    std::printf("%-52s %10.2e  %s\n", what, val, ok ? "ok" : "FAIL");
    if (!ok) ++failures;
}

int main() {
    // ---- associated Legendre functions (Condon-Shortley phase) ----
    {
        double e = 0;
        for (double th : {0.0, 0.3, 1.1, 2.0, std::numbers::pi}) {
            const double mu = std::cos(th), s = std::sin(th);
            e = std::max(e, std::abs(legendre(2, 0, th).P - (3 * mu * mu - 1) / 2));
            e = std::max(e, std::abs(legendre(2, 1, th).P - (-3 * mu * s)));
            e = std::max(e, std::abs(legendre(2, 2, th).P - 3 * s * s));
            e = std::max(e, std::abs(legendre(3, 1, th).P - (-1.5 * (5 * mu * mu - 1) * s)));
            e = std::max(e, std::abs(legendre(3, 1, th).Ps - (-1.5 * (5 * mu * mu - 1))));
            e = std::max(e, std::abs(legendre(4, 3, th).P - (-105 * mu * s * s * s)));
        }
        check(e < 1e-12, "Legendre: closed forms", e);
        double ed = 0, h = 1e-6;
        for (int l = 2; l <= 6; ++l)
            for (int m = 0; m <= l; ++m)
                for (double th : {0.2, 0.9, 1.7, 2.8}) {
                    const double fd =
                        (legendre(l, m, th + h).P - legendre(l, m, th - h).P) / (2 * h);
                    ed = std::max(ed, std::abs(legendre(l, m, th).dP - fd) /
                                          (1 + std::abs(fd)));
                    if (m > 0)
                        ed = std::max(ed, std::abs(legendre(l, m, th).Ps * std::sin(th) -
                                                   legendre(l, m, th).P));
                }
        check(ed < 1e-7, "Legendre: dP/dtheta and P/sin(theta)", ed);
    }

    // ---- eigenfunctions ----
    struct Case { double Oh; int l, m; };
    for (Case cs : {Case{std::sqrt(1.60997e-6), 2, 0}, Case{0.0316, 3, 2}, Case{0.316, 4, 1}}) {
        const ModeResult r = solve_mode(cs.l, cs.Oh, 48, std::nullopt, 1e-10);
        const Eigenfunction ef(r, cs.m);
        const cplx lam = r.lambda;
        char buf[96];

        // kinematic condition u(1) = lambda d, d = 1
        const double ek = std::abs(ef.radial(1.0).u - lam) / std::abs(lam);
        std::snprintf(buf, sizeof buf, "Oh %.3g l %d: kinematic u(1) = lambda", cs.Oh, cs.l);
        check(ek < 1e-12, buf, ek);

        // tangential stress v' - v + u = 0 at r = 1 (one-sided finite difference)
        const double h = 1e-5;
        const cplx v1 = (3.0 * ef.radial(1.0).v - 4.0 * ef.radial(1 - h).v + ef.radial(1 - 2 * h).v) / (2 * h);
        const cplx a1 = v1 - ef.radial(1.0).v + ef.radial(1.0).u;
        const double et = std::abs(a1) / std::abs(lam);
        std::snprintf(buf, sizeof buf, "Oh %.3g l %d: tangential stress at r = 1", cs.Oh, cs.l);
        check(et < 1e-6, buf, et);

        // div u = 0 at an interior point (central differences in r, theta, phi)
        const double rr = 0.7, th = 1.0, ph = 0.4, d = 1e-5;
        auto F = [&](double r_, double t_, double p_) { return ef.field(r_, t_, p_); };
        const cplx div =
            (std::pow(rr + d, 2) * F(rr + d, th, ph).ur - std::pow(rr - d, 2) * F(rr - d, th, ph).ur) /
                (2 * d * rr * rr) +
            (std::sin(th + d) * F(rr, th + d, ph).ut - std::sin(th - d) * F(rr, th - d, ph).ut) /
                (2 * d * rr * std::sin(th)) +
            (F(rr, th, ph + d).uf - F(rr, th, ph - d).uf) / (2 * d * rr * std::sin(th));
        const double ed = std::abs(div) / std::abs(lam);
        std::snprintf(buf, sizeof buf, "Oh %.3g l %d m %d: div u = 0", cs.Oh, cs.l, cs.m);
        check(ed < 1e-6, buf, ed);

        // small Oh: potential flow u = lambda r^(l-1), v = lambda r^(l-1)/l away from the surface
        if (cs.Oh < 1e-2) {
            const double r0 = 0.5;
            const auto a = ef.radial(r0);
            const cplx up = lam * std::pow(r0, cs.l - 1);
            const double ep = std::max(std::abs(a.u - up), std::abs(a.v - up / double(cs.l))) /
                              std::abs(up);
            check(ep < 1e-2, "small Oh: potential flow at r = 0.5", ep);
        }
    }
    std::printf(failures ? "FAIL\n" : "ok\n");
    return failures ? 1 : 0;
}
