/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
