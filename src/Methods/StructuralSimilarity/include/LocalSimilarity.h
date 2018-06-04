//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_LOCALSIMILARITY_H
#define AMOLQCPP_LOCALSIMILARITY_H

#include <complex>
#include "NeighborhoodExpander.h"

namespace LocalSimilarity {

    double unnormalizedLocalSimilarity(const Environment &e1,
                                       const Environment &e2);

    double localSimilarity(const Environment& e1,
                           const Environment& e2,
                           double zeta = ExpansionSettings::zeta);

    double unnormalizedLocalSelfSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions);
    double unnormalizedLocalSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions1,
                                       const TypeSpecificNeighborhoodsAtOneCenter& expansions2);

    double localSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions1,
                           const TypeSpecificNeighborhoodsAtOneCenter& expansions2,
                           double zeta = ExpansionSettings::zeta);


    namespace {
        double generic(const TypeSpecificNeighborhoodsAtOneCenter &expansions);
        double generic(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                       const TypeSpecificNeighborhoodsAtOneCenter &expansions2);

        double chemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions);
        double chemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                        const TypeSpecificNeighborhoodsAtOneCenter &expansions2);

        double kroneckerDelta(int typeA, int typeB);
        double alchemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                          const TypeSpecificNeighborhoodsAtOneCenter &expansions2);
    }

    //std::complex<double> distance(const PowerSpectrum &a, const PowerSpectrum &b)
};

#endif //AMOLQCPP_LOCALSIMILARITY_H
