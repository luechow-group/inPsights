//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_ENVIRONMENTEXPANDER_H
#define AMOLQCPP_ENVIRONMENTEXPANDER_H

#include "RadialGaussianBasis.h"
#include "NeighborhoodExpansion.h"
#include "Environment.h"
#include <Type.h>

//using TypemappedNeighborhoods = std::map<int,NeighborhoodExpansion>; // generic would be a single on
////using ParticularNeighborhoods = std::vector<TypemappedNeighborhood>; // contains all expansions of one particle center
//// map containing all centers of the molecule (according to the molecules atom + electron order)
//using MolecularNeighborhoods = std::map<NumberedType<int>,TypemappedNeighborhoods>;

using TypeSpecificNeighborhoodsAtOneCenter = std::map<int, NeighborhoodExpansion>; // expansion coeffs related to an int type

using MolecularCenters = std::map<NumberedType<int>,TypeSpecificNeighborhoodsAtOneCenter>;

class NeighborhoodExpander{
public:
    explicit NeighborhoodExpander();

    std::complex<double>
    coefficient(unsigned n, unsigned l, int m, const SphericalCoordinates &coords, double weight, double weightScale,
                    double neighborSigma = ExpansionSettings::Radial::sigmaAtom) const;

    NeighborhoodExpansion expandEnvironment(const Environment& e, int expansionTypeId = 0) const;

    TypeSpecificNeighborhoodsAtOneCenter computeExpansions(const Environment &e);

    MolecularCenters computeExpansions(MolecularGeometry molecule);

private:
    RadialGaussianBasis radialGaussianBasis_;
};

#endif //AMOLQCPP_ENVIRONMENTEXPANDER_H
