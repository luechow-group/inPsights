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
