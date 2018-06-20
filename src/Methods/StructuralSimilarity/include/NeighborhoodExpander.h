//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_ENVIRONMENTEXPANDER_H
#define AMOLQCPP_ENVIRONMENTEXPANDER_H

#include "RadialGaussianBasis.h"
#include "NeighborhoodExpansion.h"
#include "Environment.h"

using TypeSpecificNeighborhoodsAtOneCenter = std::map<int, NeighborhoodExpansion>; // expansion coeffs related to an int type

using MolecularCenters = std::map<NumberedType<int>,TypeSpecificNeighborhoodsAtOneCenter>;

class NeighborhoodExpander{
public:
    explicit NeighborhoodExpander();

    NeighborhoodExpansion expandEnvironment(const Environment& e, int expansionTypeId = 0) const;

    TypeSpecificNeighborhoodsAtOneCenter computeParticularExpansions(const Environment &e);

    MolecularCenters computeMolecularExpansions(MolecularGeometry molecule);

private:
    RadialGaussianBasis radialGaussianBasis_;
};

#endif //AMOLQCPP_ENVIRONMENTEXPANDER_H
