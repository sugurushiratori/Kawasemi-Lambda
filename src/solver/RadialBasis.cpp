#include "solver/RadialBasis.h"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace lambda {

Eigen::VectorXd RadialBasis::points(int Nr) {
    // r_j = cos(pi/2 (1 - j/(Nr-1))): Gauss-Lobatto points of [-1, 1] restricted to r > 0;
    // the point r = 0 is replaced by the midpoint of its neighbourhood (j = 0).
    Eigen::VectorXd r(Nr);
    for (int j = 0; j < Nr; ++j) {
        double s = static_cast<double>(j) / (Nr - 1);
        if (j == 0) s = 1.0 / (2.0 * (Nr - 1));
        r[j] = std::cos(std::numbers::pi / 2.0 * (1.0 - s));
    }
    r[Nr - 1] = 1.0;
    return r;
}

void RadialBasis::evaluate(int Nr, int parity, double x, Eigen::VectorXd& T, Eigen::VectorXd& T1,
                           Eigen::VectorXd& T2) {
    const int K = 2 * Nr;  // degrees 0 .. 2Nr-1
    Eigen::VectorXd t(K), t1(K), t2(K);
    // T_{k+1} = 2x T_k - T_{k-1}, differentiated once and twice (stable at x = 1)
    t[0] = 1.0;  t1[0] = 0.0;  t2[0] = 0.0;
    t[1] = x;    t1[1] = 1.0;  t2[1] = 0.0;
    for (int k = 1; k + 1 < K; ++k) {
        t[k + 1] = 2.0 * x * t[k] - t[k - 1];
        t1[k + 1] = 2.0 * t[k] + 2.0 * x * t1[k] - t1[k - 1];
        t2[k + 1] = 4.0 * t1[k] + 2.0 * x * t2[k] - t2[k - 1];
    }
    T.resize(Nr);
    T1.resize(Nr);
    T2.resize(Nr);
    for (int n = 0; n < Nr; ++n) {
        T[n] = t[2 * n + parity];
        T1[n] = t1[2 * n + parity];
        T2[n] = t2[2 * n + parity];
    }
}

RadialBasis::RadialBasis(int Nr, int parity) : Nr_(Nr), parity_(parity) {
    if (Nr < 4) throw std::invalid_argument("Nr must be >= 4");
    const Eigen::VectorXd r = points(Nr);
    T_.resize(Nr, Nr);
    T1_.resize(Nr, Nr);
    T2_.resize(Nr, Nr);
    Eigen::VectorXd t, t1, t2;
    for (int j = 0; j < Nr; ++j) {
        evaluate(Nr, parity, r[j], t, t1, t2);
        T_.col(j) = t;
        T1_.col(j) = t1;
        T2_.col(j) = t2;
    }
}

}  // namespace lambda
