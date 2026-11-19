// Spatial form of an eigenmode.
//
// With the surface displacement d = 1 (eigenvector normalised in the solver) and
// Y(theta, phi) = P_l^m(cos theta) e^{i m phi} / max|P_l^m|  (so that max|Y| = 1):
//
//   surface   r = 1 + Re[ d Y e^{lambda t} ]
//   u_r     = u(r) Y,   u_theta = v(r) dY/dtheta,   u_phi = v(r) (1/sin theta) dY/dphi
//   p       = p(r) Y
//
// All quantities are complex amplitudes; the physical field at time t is
// Re[ q e^{lambda t} ] = e^{-Gamma t} (Re q cos(omega t) - Im q sin(omega t)).
// Units: length R, velocity nu/R, pressure rho nu^2/R^2 (viscous units).
#pragma once

#include <complex>

#include "solver/Solver.h"

namespace lambda {

// Associated Legendre function P_l^m(cos theta) (Condon-Shortley phase) with
// dP/dtheta and P/sin(theta), all finite at the poles.
struct Legendre {
    double P = 0.0, dP = 0.0, Ps = 0.0;
};
Legendre legendre(int l, int m, double theta);

struct RadialAmplitude {
    cplx u, v, p;
};

struct FieldAmplitude {
    cplx ur, ut, uf, p;
};

class Eigenfunction {
public:
    Eigenfunction(const ModeResult& mode, int m);

    int l() const { return l_; }
    int m() const { return m_; }
    cplx lambda() const { return lambda_; }

    RadialAmplitude radial(double r) const;
    FieldAmplitude field(double r, double theta, double phi) const;
    cplx surface(double theta, double phi) const;  // d Y

private:
    int l_, m_, Nr_;
    cplx lambda_;
    Eigen::VectorXcd c_;
    double ynorm_ = 1.0;  // 1 / max|P_l^m|
};

}  // namespace lambda
