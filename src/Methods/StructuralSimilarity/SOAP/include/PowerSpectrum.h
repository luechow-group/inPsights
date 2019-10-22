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

#ifndef INPSIGHTS_POWERSPECTRUM_H
#define INPSIGHTS_POWERSPECTRUM_H

#include <Eigen/Core>
#include "SOAPSettings.h"

namespace SOAP {
    class NeighborhoodExpansion;

    namespace PowerSpectrum {
        //function to calculate p_ab(X_i)
        Eigen::VectorXcd partialPowerSpectrum(const NeighborhoodExpansion &n1a,
                                              const NeighborhoodExpansion &n1b);

        Eigen::VectorXcd partialPowerSpectrum(const NeighborhoodExpansion &n); //TODO CAN THIS SAFE COMPUTATIONAL TIME?

        // PARTICLE-TYPE SPECIFIC
        std::complex<double> powerSpectrumCoefficient(const NeighborhoodExpansion &speciesA,/*replace by Type typeA*/
                                                      const NeighborhoodExpansion &speciesB,/*replace by Type typeB*/
                                                      unsigned n1, unsigned n2, unsigned l);

        // GENERIC
        std::complex<double>
        powerSpectrumCoefficient(const NeighborhoodExpansion &generic,/*replace GenericType generic*/
                                 unsigned n1, unsigned n2, unsigned l);

    }
}

#endif //INPSIGHTS_POWERSPECTRUM_H
