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
                           unsigned zeta = ExpansionSettings::zeta);

    //std::complex<double> distance(const PowerSpectrum &a, const PowerSpectrum &b)

    double unnormalizedLocalSimilarity(const std::map<int, NeighborhoodExpansion> &expansions1,
                                       const std::map<int, NeighborhoodExpansion> &expansions2);

    double localSimilarity(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                           const TypeSpecificNeighborhoodsAtOneCenter &expansions2,
                           double zeta = ExpansionSettings::zeta);
};

#endif //AMOLQCPP_LOCALSIMILARITY_H
