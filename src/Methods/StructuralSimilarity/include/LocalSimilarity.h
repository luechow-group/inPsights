//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_LOCALSIMILARITY_H
#define AMOLQCPP_LOCALSIMILARITY_H

#include <complex>
#include "NeighborhoodExpander.h"

namespace LocalSimilarity {

    double unnormalizedLocalSimialrity(const Environment& e1,
                                       const Environment& e2);
    double localSimilarity(const Environment& e1,
                           const Environment& e2,
                           double zeta = ExpansionSettings::zeta);

    double unnormalizedLocalSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions1,
                                       const TypeSpecificNeighborhoodsAtOneCenter& expansions2);

    double unnormalizedLocalSelfSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions);


    double localSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions1,
                           const TypeSpecificNeighborhoodsAtOneCenter& expansions2,
                           double zeta = ExpansionSettings::zeta);

    //std::complex<double> distance(const PowerSpectrum &a, const PowerSpectrum &b)
};

#endif //AMOLQCPP_LOCALSIMILARITY_H
