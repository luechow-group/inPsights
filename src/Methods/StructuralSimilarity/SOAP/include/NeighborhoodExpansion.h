// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_NEIGHBORHOODEXPANSION_H
#define INPSIGHTS_NEIGHBORHOODEXPANSION_H

#include <Eigen/Core>

// contains the expansions of all particles with identical type
// that were within a given cutoff radius around a single center

namespace SOAP {
    class NeighborhoodExpansion {
    public:
        // Preallocates the coefficient matrix
        NeighborhoodExpansion();

        std::complex<double> getCoefficient(unsigned n, unsigned l, int m) const;

        Eigen::Ref<const Eigen::VectorXcd> getCoefficients_nl(unsigned n, unsigned l) const;

        void storeCoefficient(unsigned n, unsigned l, int m, const std::complex<double> &coefficient);

        Eigen::VectorXcd asEigenVector() const;

        void operator*=(double weight);

        friend std::ostream &operator<<(std::ostream &os, const NeighborhoodExpansion &ne);

        unsigned angularSubEntityLength(unsigned l) const;

        unsigned angularEntityLength(int l) const;

    private:
        unsigned angularEntityLength_, entityLength_;
        Eigen::VectorXcd coefficients_;
    };
}

#endif //INPSIGHTS_NEIGHBORHOODEXPANSION_H
