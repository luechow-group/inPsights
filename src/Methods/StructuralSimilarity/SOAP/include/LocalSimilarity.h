//
// Created by Michael Heuer on 09.05.18.
//

#ifndef INPSIGHTS_LOCALSIMILARITY_H
#define INPSIGHTS_LOCALSIMILARITY_H

#include <complex>
#include "NeighborhoodExpander.h"
#include "ExpansionSettings.h"

namespace LocalSimilarity {

    double unnormalizedKernel(const Environment &e1,
                              const Environment &e2);

    double unnormalizedSelfKernel(const Environment &e);

    double kernel(const Environment &e1,
                  const Environment &e2,
                  double zeta = SOAPExpansion::settings.zeta.get());

    double unnormalizedSelfKernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions);
    double unnormalizedKernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                              const TypeSpecificNeighborhoodsAtOneCenter &expansions2);

    double kernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                  const TypeSpecificNeighborhoodsAtOneCenter &expansions2,
                  double zeta = SOAPExpansion::settings.zeta.get());

    double kernelDistance(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                          const TypeSpecificNeighborhoodsAtOneCenter &expansions2,
                          double zeta = SOAPExpansion::settings.zeta.get());

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
};

#endif //INPSIGHTS_LOCALSIMILARITY_H
