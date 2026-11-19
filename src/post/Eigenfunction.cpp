#include "post/Eigenfunction.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "solver/RadialBasis.h"

namespace lambda {

Legendre legendre(int l, int m, double theta) {
    if (m < 0 || m > l) throw std::invalid_argument("need 0 <= m <= l");
    const double mu = std::cos(theta), s = std::sin(theta);
    // D[k][n] = d^k P_n / dmu^k, from (n+1) P_{n+1} = (2n+1) mu P_n - n P_{n-1}
    // differentiated k times.
    const int K = m + 1;
    std::vector<std::vector<double>> D(K + 1, std::vector<double>(l + 2, 0.0));
    for (int k = 0; k <= K; ++k) {
        double pm1 = 0.0, p0 = (k == 0) ? 1.0 : 0.0;
        D[k][0] = p0;
        for (int n = 0; n < l + 1; ++n) {
            const double dk1 = (k > 0) ? D[k - 1][n] : 0.0;
            const double p1 = ((2 * n + 1) * (mu * p0 + k * dk1) - n * pm1) / (n + 1);
            D[k][n + 1] = p1;
            pm1 = p0;
            p0 = p1;
        }
    }
    const double sign = (m % 2) ? -1.0 : 1.0;
    Legendre L;
    if (m == 0) {
        L.P = D[0][l];
        L.dP = -s * D[1][l];
        L.Ps = 0.0;  // only used multiplied by m
    } else {
        const double sm1 = std::pow(s, m - 1);
        L.P = sign * sm1 * s * D[m][l];
        L.Ps = sign * sm1 * D[m][l];
        L.dP = sign * (m * sm1 * mu * D[m][l] - sm1 * s * s * D[m + 1][l]);
    }
    return L;
}

Eigenfunction::Eigenfunction(const ModeResult& mode, int m)
    : l_(mode.l), m_(m), Nr_(mode.Nr), lambda_(mode.lambda), c_(mode.coeffs) {
    if (c_.size() != 3 * Nr_ + 1) throw std::runtime_error("no eigenvector stored for this mode");
    if (m < 0 || m > l_) throw std::invalid_argument("eigenfunction: need 0 <= m <= l");
    double pmax = 0.0;
    for (int i = 0; i <= 4000; ++i)
        pmax = std::max(pmax, std::abs(legendre(l_, m_, std::numbers::pi * i / 4000.0).P));
    ynorm_ = 1.0 / pmax;
}

RadialAmplitude Eigenfunction::radial(double r) const {
    Eigen::VectorXd V, V1, V2, S, S1, S2;
    RadialBasis::evaluate(Nr_, (l_ + 1) % 2, r, V, V1, V2);
    RadialBasis::evaluate(Nr_, l_ % 2, r, S, S1, S2);
    RadialAmplitude a;
    a.u = a.v = a.p = 0.0;
    for (int n = 0; n < Nr_; ++n) {
        a.u += V[n] * c_[n];
        a.v += V[n] * c_[Nr_ + n];
        a.p += S[n] * c_[2 * Nr_ + n];
    }
    return a;
}

FieldAmplitude Eigenfunction::field(double r, double theta, double phi) const {
    const RadialAmplitude a = radial(r);
    const Legendre L = legendre(l_, m_, theta);
    const cplx e = std::polar(ynorm_, m_ * phi);
    FieldAmplitude f;
    f.ur = a.u * L.P * e;
    f.ut = a.v * L.dP * e;
    f.uf = a.v * cplx(0.0, m_) * L.Ps * e;
    f.p = a.p * L.P * e;
    return f;
}

cplx Eigenfunction::surface(double theta, double phi) const {
    const Legendre L = legendre(l_, m_, theta);
    return c_[3 * Nr_] * L.P * std::polar(ynorm_, m_ * phi);
}

}  // namespace lambda
