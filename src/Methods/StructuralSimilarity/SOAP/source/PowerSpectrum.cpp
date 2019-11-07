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

#include "PowerSpectrum.h"
#include "NeighborhoodExpansion.h"

using namespace SOAP;

//function to calculate p_ab(X_i)
Eigen::VectorXcd PowerSpectrum::partialPowerSpectrum(const NeighborhoodExpansion& n1a,
                                                     const NeighborhoodExpansion& n1b) {
    const auto nmax = Radial::settings.nmax();
    const auto lmax = Angular::settings.lmax();

    unsigned angularEntityLength = 2 * lmax + 1;
    unsigned entityLength = nmax * nmax * angularEntityLength;
    Eigen::VectorXcd expansionCoefficients = Eigen::VectorXcd::Zero(entityLength);

    for (unsigned n1 = 1; n1 <= nmax; ++n1) {
        unsigned n1BlockStartIdx = (n1 - 1) * nmax * angularEntityLength;

        for (unsigned n2 = 1; n2 <= nmax; ++n2) {
            unsigned n1n2BlockStartIdx = n1BlockStartIdx + (n2 - 1) * angularEntityLength;

            for (unsigned l = 0; l <= lmax; ++l) {
                //TODO NORM HERE?
                //expansionCoefficients[n1n2BlockStartIdx + l] = std::norm(powerSpectrumCoefficient(n1a, n1b,n1, n2, l));
                expansionCoefficients[n1n2BlockStartIdx + l] = powerSpectrumCoefficient(n1a, n1b,n1, n2, l);
            }
        }
    }
    return expansionCoefficients;
}

std::complex<double> PowerSpectrum::powerSpectrumCoefficient(const NeighborhoodExpansion& generic,
                                                             unsigned n1, unsigned n2, unsigned l) {
    return powerSpectrumCoefficient(generic,generic,n1,n2,l); // TODO CAN WE SAVE EFFORT HERE?
}

std::complex<double> PowerSpectrum::powerSpectrumCoefficient(const NeighborhoodExpansion& speciesA,
                                                             const NeighborhoodExpansion& speciesB,
                                                             unsigned n1, unsigned n2, unsigned l ) {
    Radial::checkBounds(n1);
    Radial::checkBounds(n2);
    Angular::checkBounds(l);

    std::complex<double> sum = (speciesA.getCoefficients_nl(n1, l).array().conjugate()
                  * speciesB.getCoefficients_nl(n2, l).array()).sum();
    return M_PI * sqrt(8./(2.*l+1)) * sum; // in the soapxx code, this prefactor can be replaced by 1 (using a bool)
}

