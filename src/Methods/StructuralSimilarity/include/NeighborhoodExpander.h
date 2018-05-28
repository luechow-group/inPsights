//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_ENVIRONMENTEXPANDER_H
#define AMOLQCPP_ENVIRONMENTEXPANDER_H

#include "RadialGaussianBasis.h"
#include "NeighborhoodExpansion.h"
#include "Environment.h"
#include <Type.h>

using TypeSpecificExpansionsVector = std::map<int, NeighborhoodExpansion>; // expansion coeffs related to an int type

using AllCentersSet = std::vector<TypeSpecificExpansionsVector>;

class NeighborhoodExpander{
public:
    explicit NeighborhoodExpander();

    std::complex<double>
    coefficient(unsigned n, unsigned l, int m, const SphericalCoordinates &coords, double weight, double weightScale,
                    double neighborSigma = ExpansionSettings::Radial::sigmaAtom) const;

    NeighborhoodExpansion expandEnvironment(const Environment& e, int expansionTypeId = 0) const;

    TypeSpecificExpansionsVector computeExpansions(const Environment &e);

    AllCentersSet computeExpansions(MolecularGeometry molecule, GeneralStorageType type);//OTHER TYPE!! RENAME

private:
    RadialGaussianBasis radialGaussianBasis_;
};

#endif //AMOLQCPP_ENVIRONMENTEXPANDER_H
