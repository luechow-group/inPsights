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
    explicit CoefficientsVector(unsigned numberOfParticles,
                                const ExpansionSettings& settings = ExpansionSettings::defaults())
            : s_(settings),
              entityLength_(unsigned( s_.radial.nmax * (pow(s_.angular.lmax,2) + 2*s_.angular.lmax + 1) )),
              coefficients_(Eigen::VectorXcd::Zero(numberOfParticles * entityLength_))
    {
        // increment counter for each entity
        for (int j = 0; j < numberOfParticles; ++j) incrementNumberOfEntities();
    }

    explicit CoefficientsVector(const ExpansionSettings& settings = ExpansionSettings::defaults())
            : CoefficientsVector(1, settings)
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

    void operator*=(double weight){
        coefficients_ *= weight;
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
