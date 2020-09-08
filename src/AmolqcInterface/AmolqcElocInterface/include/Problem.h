// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef PROBLEM_H
#define PROBLEM_H

#include <array>
#include <vector>
#include <Eigen/Core>

//#include "Meta.h"

namespace cppoptlib {

template<typename Scalar_, int Dim_ = Eigen::Dynamic>
class Problem {
public:
    static const int Dim = Dim_;
    typedef Scalar_ Scalar;
    using TVector   = Eigen::Matrix<Scalar, Dim, 1>;
    using THessian  = Eigen::Matrix<Scalar, Dim, Dim>;
    //using TCriteria = Criteria<Scalar>;
    using TIndex = typename TVector::Index;
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    bool gradientResetQ{false};

public:
    Problem() = default;

    virtual ~Problem() = default;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"


#pragma GCC diagnostic pop

    virtual Scalar value(const TVector &x) = 0;

    Scalar operator()(const TVector &x) {
        return value(x);
    }

    virtual void gradient(const TVector &x, TVector &grad) {
        finiteGradient(x, grad);
    }

    Scalar calculateFirstDerivativeDeltaNash(Scalar xi) {
        auto sqrtEps = std::sqrt(std::numeric_limits<Scalar>::epsilon());
        //from John C. Nash Compact Numerical Methods for Computers 2ed, p. 219
        return sqrtEps * (std::abs(xi) + sqrtEps) *
               1E2; // *1E2 as a safety margin - rather loose digits than having a too small h
    }

    Scalar calculateSecondDerivativeDeltaNash(Scalar gradi) {
        auto cbrtEps = std::cbrt(std::numeric_limits<Scalar>::epsilon());
        // from https://mathoverflow.net/questions/28463/optimum-small-number-for-numerical-differentiation
        return cbrtEps * (std::abs(gradi) + cbrtEps) *
               1E3; // *1E3 as a safety margin - rather loose digits than having a too small h
    }

    virtual void hessian(const TVector &x, THessian &hessian) {
        finiteHessian(x, hessian);
    }

    virtual bool checkGradient(const TVector &x, int accuracy = 3) {
        const TIndex D = x.rows();
        TVector actual_grad(D);
        TVector expected_grad(D);
        gradient(x, actual_grad);
        finiteGradient(x, expected_grad, accuracy);
        for (TIndex d = 0; d < D; ++d) {
            Scalar scale = std::max(static_cast<Scalar>(std::max(fabs(actual_grad[d]), fabs(expected_grad[d]))),
                                    Scalar(1.));
            if (fabs(actual_grad[d] - expected_grad[d]) > 1e-2 * scale)
                return false;
        }
        return true;

    }

    virtual bool checkHessian(const TVector &x, int accuracy = 3) {
        const TIndex D = x.rows();

        THessian actual_hessian = THessian::Zero(D, D);
        THessian expected_hessian = THessian::Zero(D, D);
        hessian(x, actual_hessian);
        finiteHessian(x, expected_hessian, accuracy);
        for (TIndex d = 0; d < D; ++d) {
            for (TIndex e = 0; e < D; ++e) {
                Scalar scale = std::max(
                        static_cast<Scalar>(std::max(fabs(actual_hessian(d, e)), fabs(expected_hessian(d, e)))),
                        Scalar(1.));
                if (fabs(actual_hessian(d, e) - expected_hessian(d, e)) > 1e-1 * scale)
                    return false;
            }
        }
        return true;
    }

    void finiteGradient(const TVector &x, TVector &grad, int accuracy = 0) {
        // accuracy can be 0, 1, 2, 3
        //const Scalar eps = 2.2204e-6; // Original Implementation
        const TIndex D = x.rows();
        static const std::array<std::vector<Scalar>, 4> coeff =
                {{{1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3}}};
        static const std::array<std::vector<Scalar>, 4> coeff2 =
                {{{1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4}}};
        static const std::array<Scalar, 4> dd = {2, 12, 60, 840};

        grad.resize(D);
        TVector &xx = const_cast<TVector &>(x);

        const int innerSteps = 2 * (accuracy + 1);
        //const Scalar ddVal = dd[accuracy]*eps; // Original Implementation

        for (TIndex d = 0; d < D; d++) {
            grad[d] = 0;
            Scalar eps = calculateFirstDerivativeDeltaNash(x[d]); // Modified
            Scalar ddVal = dd[accuracy] * eps; // Modified
            for (int s = 0; s < innerSteps; ++s) {
                Scalar tmp = xx[d];
                xx[d] += coeff2[accuracy][s] * eps;
                grad[d] += coeff[accuracy][s] * value(xx);
                xx[d] = tmp;
            }
            grad[d] /= ddVal;
        }
    }

    void riddersFiniteGradient(const TVector &x, TVector &grad, int accuracy = 0) { //Numerical Recipes
        const TIndex D = x.rows();
        Scalar err;
        const int ntab = 10; //Sets maximum size of tableau.
        const Scalar con = 1.4, con2 = (con * con); //Stepsize decreased by CON at each iteration.
        const Scalar big = std::numeric_limits<Scalar>::max();
        const Scalar safe = 2.0; //Return when error is SAFE worse than the
        int i, j; //best so far.
        Scalar errt = big, fac = 0, hh = 0, ans = 0;
        Eigen::Matrix<Scalar, ntab, ntab> a;

        TVector xiphh, ximhh;

        for (TIndex d = 0; d < x.rows(); d++) {
            Scalar h = calculateFirstDerivativeDeltaNash(x[d]);

            assert(h != 0.0 && "h must be nonzero.");
            hh = h;
            xiphh = x;
            xiphh(d) += hh;
            ximhh = x;
            ximhh(d) -= hh;
            a(0, 0) = (value(xiphh) - value(ximhh)) / (2.0 * hh);
            err = big;
            for (i = 1; i < ntab; i++) {
                //Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
                hh /= con;
                xiphh = x;
                xiphh(d) += hh;
                ximhh = x;
                ximhh(d) -= hh;
                a(0, i) = (value(xiphh) - value(ximhh)) / (2.0 * hh); //Try new, smaller stepsize.
                fac = con2;
                for (j = 1; j <= i; j++) { //Compute extrapolations of various orders, requiring

                    a(j, i) = (a(j - 1, i) * fac - a(j - 1, i - 1)) / (fac - 1.0); //no new function evaluations.
                    fac = con2 * fac;
                    errt = std::max({std::abs(a(j, i) - a(j - 1, i)), std::abs(a(j, i) - a(j - 1, i - 1))});
                    //The error strategy is to compare each new extrapolation to one order lower, both at the present stepsize and the previous one.

                    if (errt <= err) {
                        //If error is decreased, save the improved answer.
                        err = errt;
                        ans = a(j, i);
                    }
                    //std::cout <<"errt " << errt << std::endl;
                }

                if (abs(a(i, i) - a(i - 1, i - 1)) >= safe * err) break;
                //If higher order is worse by a significant factor SAFE, then quit early.
            }
            xiphh = x;
            ximhh = x;
            //std::cout <<"err " << err << std::endl;
            grad(d) = ans;
        }
    }

    void semifiniteHessian(const TVector &x, THessian &hessian, int accuracy = 0,
                           Scalar hmax = std::numeric_limits<Scalar>::max()) {
        const TIndex D = x.rows();
        static const std::array<std::vector<Scalar>, 4> coeff =
                {{{1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3}}};
        static const std::array<std::vector<Scalar>, 4> coeff2 =
                {{{1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4}}};
        static const std::array<Scalar, 4> dd = {2, 12, 60, 840};

        auto &xx = const_cast<TVector &>(x);
        const int innerSteps = 2 * (accuracy + 1);

        TVector grad(D);
        TVector gradxx(D);
        gradient(x, grad);// for optimal h calculation

        hessian.resize(D, D);

        for (TIndex i = 0; i < D; i++) {
            Scalar h = std::min(calculateSecondDerivativeDeltaNash(grad(i)),hmax);
            for (TIndex j = 0; j < D; j++) {
                hessian(i, j) = 0;

                Scalar ddVal = dd[accuracy] * h;

                for (int s = 0; s < innerSteps; ++s) {
                    Scalar tmp = xx[j];
                    xx[j] += coeff2[accuracy][s] * h;

                    gradient(xx, gradxx);
                    hessian(i, j) += coeff[accuracy][s] * gradxx(i);
                    xx[j] = tmp;
                }
                hessian(i, j) /= ddVal;
            }
        }

        // symmetrize matrix by averaging the off-diagonal elements
        for (TIndex i = 0; i < D; i++) {
            for (TIndex j = i + 1; j < D; j++) {
                hessian(i, j) = (hessian(i, j) + hessian(j, i)) / 2.0;
                hessian(j, i) = hessian(i, j);
            }
        }
    }


    void semifiniteHessian(const TVector &x, THessian &hessian, const std::vector<unsigned long> &segmentsToCalculate,
                           unsigned long indicesPerSegment = 3, int accuracy = 0) {
        const TIndex D = x.rows();
        static const std::array<std::vector<Scalar>, 4> coeff =
                {{{1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3}}};
        static const std::array<std::vector<Scalar>, 4> coeff2 =
                {{{1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4}}};
        static const std::array<Scalar, 4> dd = {2, 12, 60, 840};

        auto &xx = const_cast<TVector &>(x);
        const int innerSteps = 2 * (accuracy + 1);

        TVector grad(D);
        TVector gradxx(D);

        gradient(x, grad);
        hessian = THessian::Zero(D, D);

        for (const auto si : segmentsToCalculate) {
            for (unsigned long i = si*3; i < si*3 + 3; ++i) {

                Scalar h = calculateSecondDerivativeDeltaNash(grad(i));
                for (const auto sj : segmentsToCalculate) {
                    for (unsigned long j = sj*3; j < sj*3 + 3; ++j) {
                        hessian(i, j) = 0;
                        Scalar ddVal = dd[accuracy] * h;

                        for (int s = 0; s < innerSteps; ++s) {
                            Scalar tmp = xx[j];
                            xx[j] += coeff2[accuracy][s] * h;

                            gradient(xx, gradxx);
                            hessian(i, j) += coeff[accuracy][s] * gradxx(i);
                            xx[j] = tmp;
                        }
                        hessian(i, j) /= ddVal;
                    }
                }
            }
        }
        // symmetrize matrix by averaging the off-diagonal elements
        for (TIndex i = 0; i < D; i++) {
            for (TIndex j = i + 1; j < D; j++) {
                hessian(i, j) = (hessian(i, j) + hessian(j, i)) / 2.0;
                hessian(j, i) = hessian(i, j);
            }
        }
        for (const auto si : segmentsToCalculate) {
            for (unsigned long i = si*3; i < si*3 + 3; ++i) {
                Scalar h = calculateSecondDerivativeDeltaNash(grad(i));
                for (const auto sj : segmentsToCalculate) {
                    for (unsigned long j = sj*3; j < sj*3 + 3; ++j) {

                        hessian(i, j) = (hessian(i, j) + hessian(j, i)) / 2.0;
                        hessian(j, i) = hessian(i, j);
                    }
                }
            }
        }
    }

  void finiteHessian(const TVector &x, THessian &hessian, int accuracy = 0) {
    const Scalar eps = std::numeric_limits<Scalar>::epsilon()*10e7;

    hessian.resize(x.rows(), x.rows());
    auto& xx = const_cast<TVector&>(x);

    if(accuracy == 0) {
      for (TIndex i = 0; i < x.rows(); i++) {
        for (TIndex j = 0; j < x.rows(); j++) {
          Scalar tmpi = xx[i];
          Scalar tmpj = xx[j];

          Scalar f4 = value(xx);
          xx[i] += eps;
          xx[j] += eps;
          Scalar f1 = value(xx);
          xx[j] -= eps;
          Scalar f2 = value(xx);
          xx[j] += eps;
          xx[i] -= eps;
          Scalar f3 = value(xx);
          hessian(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);

          xx[i] = tmpi;
          xx[j] = tmpj;
        }
      }
    } else {
      /*
        \displaystyle{{\frac{\partial^2{f}}{\partial{x}\partial{y}}}\approx
        \frac{1}{600\,h^2} \left[\begin{matrix}
          -63(f_{1,-2}+f_{2,-1}+f_{-2,1}+f_{-1,2})+\\
          63(f_{-1,-2}+f_{-2,-1}+f_{1,2}+f_{2,1})+\\
          44(f_{2,-2}+f_{-2,2}-f_{-2,-2}-f_{2,2})+\\
          74(f_{-1,-1}+f_{1,1}-f_{1,-1}-f_{-1,1})
        \end{matrix}\right] }
      */
      for (TIndex i = 0; i < x.rows(); i++) {
        for (TIndex j = 0; j < x.rows(); j++) {
          Scalar tmpi = xx[i];
          Scalar tmpj = xx[j];

          Scalar term_1 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += -2*eps;  term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += -1*eps;  term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += 1*eps;   term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += 2*eps;   term_1 += value(xx);

          Scalar term_2 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += -2*eps;  term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += -1*eps;  term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += 2*eps;   term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += 1*eps;   term_2 += value(xx);

          Scalar term_3 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += -2*eps;  term_3 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += 2*eps;   term_3 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += -2*eps;  term_3 -= value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += 2*eps;   term_3 -= value(xx);

          Scalar term_4 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += -1*eps;  term_4 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += 1*eps;   term_4 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += -1*eps;  term_4 -= value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += 1*eps;   term_4 -= value(xx);

          xx[i] = tmpi;
          xx[j] = tmpj;

          hessian(i, j) = (-63 * term_1+63 * term_2+44 * term_3+74 * term_4)/(600.0 * eps * eps);
        }
      }
    }

  }
};
}

#endif /* PROBLEM_H */
