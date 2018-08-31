//
// Created by Michael Heuer on 31.08.18.
//

#ifndef AMOLQCPP_COMPARATORS_H
#define AMOLQCPP_COMPARATORS_H

#include "Reference.h"

class ValueEuclideanDistanceComparator {
public:
    static constexpr double distThreshold = 0.01;
    static constexpr double valueThreshold = 1e-4;
    
    bool operator() (const Reference& lhs, const Reference& rhs) const;

};

#endif //AMOLQCPP_COMPARATORS_H
