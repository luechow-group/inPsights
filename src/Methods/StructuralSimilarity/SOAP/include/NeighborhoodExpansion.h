//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_NEIGHBORHOODEXPANSION_H
#define AMOLQCPP_NEIGHBORHOODEXPANSION_H

#include <Eigen/Core>

// contains the expansions of all particles with identical type
// that were within a given cutoff radius around a single center

class NeighborhoodExpansion{
public:
    // Preallocates the coefficient matrix
    NeighborhoodExpansion();

    std::complex<double> getCoefficient(unsigned n, unsigned l, int m) const;

    void storeCoefficient(unsigned n, unsigned l, int m, const std::complex<double> &coefficient);

    Eigen::VectorXcd asEigenVector() const;

    void operator*=(double weight);

    friend std::ostream& operator<<(std::ostream& os, const NeighborhoodExpansion & ne);

    unsigned angularSubEntityLength(unsigned l) const;

    unsigned angularEntityLength(int l) const;

private:
    unsigned angularEntityLength_, entityLength_;
    Eigen::VectorXcd coefficients_;
};

#endif //AMOLQCPP_NEIGHBORHOODEXPANSION_H
