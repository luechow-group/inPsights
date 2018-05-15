//
// Created by Michael Heuer on 10.05.18.
//

#ifndef AMOLQCPP_POWERSPECTRUM_H
#define AMOLQCPP_POWERSPECTRUM_H

#include <Eigen/Core>
#include "ExpansionSettings.h"
#include "ParticlePool.h"

class NeighborhoodExpansion;

namespace PowerSpectrum {
    //function to calculate p_ab(X_i)
    Eigen::VectorXd partialPowerSpectrum(const NeighborhoodExpansion& n1a,
                                         const NeighborhoodExpansion& n1b);

    // PARTICLE-TYPE SPECIFIC
    double powerSpectrumCoefficient(const NeighborhoodExpansion& speciesA,/*replace by Type typeA*/
                                    const NeighborhoodExpansion& speciesB,/*replace by Type typeB*/
                                    unsigned n1, unsigned n2, unsigned l );
    // GENERIC
    double powerSpectrumCoefficient(const NeighborhoodExpansion& generic,/*replace GenericType generic*/
                                    unsigned n1, unsigned n2, unsigned l);

};


#endif //AMOLQCPP_POWERSPECTRUM_H
