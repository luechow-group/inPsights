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

#ifndef INPSIGHTS_LOCALSIMILARITY_H
#define INPSIGHTS_LOCALSIMILARITY_H

#include <complex>
#include "NeighborhoodExpander.h"
#include "SOAPSettings.h"

namespace SOAP {
    namespace LocalSimilarity {

        double unnormalizedKernel(const Environment &e1,
                                  const Environment &e2);

        double unnormalizedSelfKernel(const Environment &e);

        double kernel(const Environment &e1,
                      const Environment &e2,
                      double zeta = General::settings.zeta());

        double unnormalizedSelfKernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions);

        double unnormalizedKernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                                  const TypeSpecificNeighborhoodsAtOneCenter &expansions2);

        double kernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                      const TypeSpecificNeighborhoodsAtOneCenter &expansions2,
                      double zeta = General::settings.zeta());

        double kernelDistance(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                              const TypeSpecificNeighborhoodsAtOneCenter &expansions2,
                              double zeta = General::settings.zeta());

        namespace internal {
            double typeAgnostic(const TypeSpecificNeighborhoodsAtOneCenter &expansions);

            double typeAgnostic(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                                const TypeSpecificNeighborhoodsAtOneCenter &expansions2);

            double chemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions);

            double chemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                            const TypeSpecificNeighborhoodsAtOneCenter &expansions2);

            double kroneckerDelta(int typeA, int typeB);

            double alchemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                              const TypeSpecificNeighborhoodsAtOneCenter &expansions2);
        }
    }
}

#endif //INPSIGHTS_LOCALSIMILARITY_H
