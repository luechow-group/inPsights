//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_COEFFICIENTSVECTOR_H
#define AMOLQCPP_COEFFICIENTSVECTOR_H

#include <Eigen/Core>

class CoefficientsVector{
public:

    // Preallocates the coefficient matrix
    CoefficientsVector(unsigned nmax, unsigned lmax);

    std::complex<double> getCoefficient(unsigned n, unsigned l, int m) const;

    void storeCoefficient(unsigned n, unsigned l, int m, const std::complex<double>& coefficient);

    void checkAngularBounds(unsigned l, int m) const;

    void checkRadialBounds(unsigned n, unsigned l) const;
private:
    unsigned nmax_, lmax_;
    Eigen::VectorXcd coefficients_;
};

#endif //AMOLQCPP_COEFFICIENTSVECTOR_H
