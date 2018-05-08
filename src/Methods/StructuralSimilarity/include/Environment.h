//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_ENVIRONMENT_H
#define AMOLQCPP_ENVIRONMENT_H

#include <Eigen/Core>
#include "ExpansionSettings.h"
#include <ParticlesVector.h>

template <typename Type>
class Environment : AbstractVector {
public:

    // Preallocates the coefficient matrix
    explicit Environment(unsigned numberOfParticles,
                                                     const ExpansionSettings& settings = ExpansionSettings::defaults())
            : AbstractVector(numberOfParticles),
              s_(settings),
              angularEntityLength_( s_.angular.lmax*s_.angular.lmax + 2*s_.angular.lmax + 1),
              entityLength_(s_.radial.nmax * angularEntityLength_),
              coefficients_(Eigen::VectorXcd::Zero(numberOfParticles * entityLength_)),
              powerSpectrum_(Eigen::VectorXd::Zero(s_.radial.nmax*s_.radial.nmax* (s_.angular.lmax+1)))
    {}

    std::complex<double> getCoefficient(unsigned i, unsigned n, unsigned l, int m) const {
    s_.checkBounds(n,l,m);
    return coefficients_( calculateIndex(i) + (n-1)*angularEntityLength_ + (l*l+l+m) );
}

    void storeCoefficient(unsigned i, unsigned n, unsigned l, int m, const std::complex<double> &coefficient) {
        s_.checkBounds(n,l,m);
        coefficients_( calculateIndex(i) + (n-1)*angularEntityLength_ + (l*l+l+m) ) = coefficient;
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

    void storeSingleNeighborExpansion(unsigned i,
                                      const Environment<Type> &singleNeighborExpansionCoefficients) {
        assert(s_ == singleNeighborExpansionCoefficients.s_
               && "The expansion settings must be identical.");
        assert(singleNeighborExpansionCoefficients.numberOfEntities() == 1
               && "CoeffcientsVector of the single particle must have the size 1.");
        coefficients_.segment( calculateIndex(i), entityLength_) = singleNeighborExpansionCoefficients.asEigenVector();
    }

    void operator*=(double weight){
        coefficients_ *= weight;
    }

    const ExpansionSettings& getSettings() const {
        return s_;
    }


    friend std::ostream& operator<<(std::ostream& os, const Environment & cv){

        auto s = cv.getSettings();

        for (unsigned i = 0; i < cv.numberOfEntities(); i++) {
            os << "i: " << i << std::endl;
            for (unsigned n = 1; n <= s.radial.nmax; ++n) {
                os << " n: "<< n << std::endl;
                for (unsigned l = 0; l <= s.angular.lmax; ++l) {
                    os << "  l: "<< l << std::endl;
                    for (int m = -int(l); m <= int(l); ++m) {
                       os << "   m: " << m << " " << cv.getCoefficient(i,n,l,m) << std::endl;
                    }
                }
            }
        }
        return os;
    }

    double calculatePowerSpectrumCoefficient(unsigned n1, unsigned n2, unsigned l) {

        double sum = 0;

        for (int m = -int(l); m < int(l); ++m) {
            sum += std::conj(getCoefficient(n1,l,m)) * getCoefficient(n2,l,m);
        }

        return M_PI*sqrt(8./(2.*l+1))* sum;
    }


private:
    ExpansionSettings s_;
    unsigned  angularEntityLength_, entityLength_;
    Eigen::VectorXcd coefficients_;
    Eigen::VectorXd powerSpectrum_;
};

#endif //AMOLQCPP_ENVIRONMENT_H
