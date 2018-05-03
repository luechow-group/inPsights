//
// Created by Michael Heuer on 03.05.18.
//

#include "CoefficientsVector.h"

CoefficientsVector::CoefficientsVector(unsigned nmax, unsigned lmax)
        : nmax_(nmax),
          lmax_(lmax),
          coefficients_(Eigen::VectorXcd::Zero( nmax_*(lmax_*lmax_+2*lmax_+1)))
{}

void CoefficientsVector::checkRadialBounds(unsigned n, unsigned l) const {
    assert( n <= nmax_ && "n must be smaller than nmax");
    assert( n >= 1 && "n must greater than or equal to 1");
    assert( l <= lmax_ && "l must be smaller than lmax");
}

void CoefficientsVector::checkAngularBounds(unsigned l, int m) const {
    assert( l <= lmax_ && "l must be smaller than lmax");
    assert( unsigned(abs(m)) <= lmax_ && "abs(m) must be smaller than lmax");
}

std::complex<double> CoefficientsVector::getCoefficient(unsigned n, unsigned l, int m) const{
    checkRadialBounds(n,l);
    checkAngularBounds(l,m);
    return coefficients_( (n-1) * (l*l+l+m) );
}

void CoefficientsVector::storeCoefficient(unsigned n, unsigned l, int m, const std::complex<double>& coefficient){
    checkRadialBounds(n,l);
    checkAngularBounds(l,m);
    coefficients_( (n-1) * (l*l+l+m) ) = coefficient;
}
