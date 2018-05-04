//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_COEFFICIENTSVECTOR_H
#define AMOLQCPP_COEFFICIENTSVECTOR_H

#include <Eigen/Core>
#include "ExpansionSettings.h"
#include <ParticlesVector.h>

template <typename Type>
class CoefficientsVector : AbstractVector {
public:

    // Preallocates the coefficient matrix
    explicit CoefficientsVector(ParticlesVector<Type> particleVector,
                                const ExpansionSettings& settings = ExpansionSettings::defaults())
            : s_(),
              entityLength_(unsigned( s_.radial.nmax * (pow(s_.angular.lmax,2) + 2*s_.angular.lmax + 1) )),
              coefficients_(Eigen::VectorXcd::Zero(particleVector.numberOfEntities() * entityLength_))
    {}

    explicit CoefficientsVector(Particle<Type> particleVector,
                                const ExpansionSettings& settings = ExpansionSettings::defaults())
            : s_(),
              entityLength_(unsigned( s_.radial.nmax * (pow(s_.angular.lmax,2) + 2*s_.angular.lmax + 1) )),
              coefficients_(Eigen::VectorXcd::Zero(entityLength_))
    {}

    std::complex<double> getCoefficient(unsigned i, unsigned n, unsigned l, int m) const {
    s_.checkBounds(n,l,m);
    return coefficients_( calculateIndex(i) + (n-1) * (l*l+l+m) );
}

    void storeCoefficient(unsigned i, unsigned n, unsigned l, int m, const std::complex<double> &coefficient) {
        s_.checkBounds(n,l,m);
        coefficients_( calculateIndex(i) + (n-1) * (l*l+l+m) ) = coefficient;
    }

    Eigen::VectorXcd asEigenVector() const {
        return coefficients_;
    }

    long calculateIndex(long i) const override {
        return AbstractVector::calculateIndex(i)*entityLength_;
    }

    //Atomic coefficients
    Eigen::VectorXcd operator[](long i) const {
        return coefficients_.segment(calculateIndex(i),entityLength_);
    }

    void storeParticleExpansion(unsigned i, const CoefficientsVector<Type> &coefficientsVector) {
        assert(s_ == coefficientsVector.s_
               && "The expansion settings must be identical.");
        assert(coefficientsVector.numberOfEntities() == 1
               && "CoeffcientsVector of the single particle must have the size 1.");
        coefficients_.segment( calculateIndex(i), entityLength_) = coefficientsVector.asEigenVector();
    }

    const ExpansionSettings& getSettings(){
        return s_;
    }

private:
    ExpansionSettings s_;
    unsigned entityLength_;
    Eigen::VectorXcd coefficients_;
};

#endif //AMOLQCPP_COEFFICIENTSVECTOR_H
