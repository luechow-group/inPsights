// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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

            double alchemicalTypeSimilarity(int typeA, int typeB);

            double alchemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                              const TypeSpecificNeighborhoodsAtOneCenter &expansions2);
        }
    }
}

#endif //INPSIGHTS_LOCALSIMILARITY_H
