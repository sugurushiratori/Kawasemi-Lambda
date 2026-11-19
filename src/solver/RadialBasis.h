// Chebyshev basis in r for the radial problem.
//
// r in (0, 1] is the right half of the Chebyshev interval [-1, 1]. Regularity
// at r = 0 fixes the parity of each radial function: velocity-type functions
// (u_r, tangential poloidal v) behave like r^(l-1), scalar-type functions (p)
// like r^l. The basis is therefore T_{2n+q}(r), n = 0..Nr-1, with the parity
// q chosen per variable.
#pragma once

#include <Eigen/Dense>

namespace lambda {

class RadialBasis {
public:
    RadialBasis(int Nr, int parity);

    int size() const { return Nr_; }
    int parity() const { return parity_; }

    // values at collocation point j of basis function n: T, dT/dr, d2T/dr2
    double T(int n, int j) const { return T_(n, j); }
    double T1(int n, int j) const { return T1_(n, j); }
    double T2(int n, int j) const { return T2_(n, j); }

    // collocation points r_j, j = 0..Nr-1, r_{Nr-1} = 1 (the surface)
    static Eigen::VectorXd points(int Nr);

    // T_{2n+parity}(r) and its first two derivatives, n = 0..Nr-1, at any r in [0, 1]
    static void evaluate(int Nr, int parity, double r, Eigen::VectorXd& T, Eigen::VectorXd& T1,
                         Eigen::VectorXd& T2);

private:
    int Nr_, parity_;
    Eigen::MatrixXd T_, T1_, T2_;
};

}  // namespace lambda
