// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
