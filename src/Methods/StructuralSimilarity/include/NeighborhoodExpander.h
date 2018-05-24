//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_ENVIRONMENTEXPANDER_H
#define AMOLQCPP_ENVIRONMENTEXPANDER_H

#include "RadialGaussianBasis.h"
#include "NeighborhoodExpansion.h"
#include "ElementType.h"
#include "Environment.h"

class NeighborhoodExpander{
public:
    explicit NeighborhoodExpander();

    std::complex<double>
    coefficient(unsigned n, unsigned l, int m, const SphericalCoordinates &coords, double weight, double weightScale,
                    double neighborSigma = ExpansionSettings::Radial::sigmaAtom) const;

    NeighborhoodExpansion expandEnvironment(const Environment& e, int expansionTypeId = 0) const;

    std::map<int, NeighborhoodExpansion> computeExpansions(const Environment &e);

private:
    RadialGaussianBasis radialGaussianBasis_;
};

#endif //AMOLQCPP_ENVIRONMENTEXPANDER_H
