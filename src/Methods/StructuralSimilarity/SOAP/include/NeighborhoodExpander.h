// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ENVIRONMENTEXPANDER_H
#define INPSIGHTS_ENVIRONMENTEXPANDER_H

#include "RadialBasis.h"
#include "NeighborhoodExpansion.h"
#include "Environment.h"

namespace SOAP {
    using TypeSpecificNeighborhoodsAtOneCenter = std::map<int, NeighborhoodExpansion>; // expansion coeffs related to an int type
    using MolecularCenters = std::map<EnumeratedType<int>, TypeSpecificNeighborhoodsAtOneCenter>;

    class NeighborhoodExpander {
    public:
        explicit NeighborhoodExpander();

        NeighborhoodExpansion expandEnvironment(const Environment &e, int expansionTypeId = 0) const;

        TypeSpecificNeighborhoodsAtOneCenter computeParticularExpansions(const Environment &e);

        MolecularCenters computeMolecularExpansions(MolecularGeometry molecule);

    private:
        RadialBasis radialGaussianBasis_;
    };
}

#endif //INPSIGHTS_ENVIRONMENTEXPANDER_H
