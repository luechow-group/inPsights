//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_NEIGHBORHOODEXPANSION_H
#define AMOLQCPP_NEIGHBORHOODEXPANSION_H

#include <Eigen/Core>

class NeighborhoodExpansion{
public:
    // Preallocates the coefficient matrix
    explicit NeighborhoodExpansion(unsigned numberOfParticles);

    std::complex<double> getCoefficient(/*unsigned i,*/ unsigned n, unsigned l, int m) const;
    //unsigned calculateIndex(unsigned i) const { return i*entityLength_; }

    void storeCoefficient(/*unsigned i,*/ unsigned n, unsigned l, int m, const std::complex<double> &coefficient);

    Eigen::VectorXcd asEigenVector() const;

    void operator*=(double weight);

    friend std::ostream& operator<<(std::ostream& os, const NeighborhoodExpansion & ne);

    unsigned getNumberOfParticles() const;

    unsigned angularSubEntityLength(unsigned l) const;
    unsigned angularEntityLength(int l) const;

private:
    unsigned numberOfParticles_;
    unsigned  /*angularSubEntityLength_,*/angularEntityLength_, entityLength_;
    Eigen::VectorXcd coefficients_;
};

#endif //AMOLQCPP_NEIGHBORHOODEXPANSION_H
