//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_LOCALSIMILARITY_H
#define AMOLQCPP_LOCALSIMILARITY_H

#include <complex>
#include "NeighborhoodExpander.h"

class Environment;

namespace LocalSimilarity {

    double unnormalizedLocalSimialrity(const Environment& e1, const Environment& e2);
    double localSimilarity(const Environment& e1, const Environment& e2);

    // static std::complex<double> distance(const PowerSpectrum &a, const PowerSpectrum &b)

    void computeGeneric(NeighborhoodExpander &neighborhoodExpander,
                        std::map<Elements::ElementType, NeighborhoodExpansion> &expansions1,
                        std::map<Elements::ElementType, NeighborhoodExpansion> &expansions2);


        double unnormalizedLocalSimialrity(const std::map<Elements::ElementType, NeighborhoodExpansion>& expansions1,
                                           const std::map<Elements::ElementType, NeighborhoodExpansion>& expansions2);

        std::map<Elements::ElementType, NeighborhoodExpansion>
        computeExpansions(const Environment &e);

};

#endif //AMOLQCPP_LOCALSIMILARITY_H
