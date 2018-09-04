//
// Created by Michael Heuer on 31.08.18.
//

#ifndef AMOLQCPP_COMPARATORS_H
#define AMOLQCPP_COMPARATORS_H

#include "Reference.h"

namespace Comparators {

    static bool identicalQ(const Reference &lhs, const Reference &rhs, double distThreshold = 0.01);
    static bool globallySimilarQ(const Reference &lhs, const Reference &rhs, double distThreshold = 0.01);

    class ValueEuclideanDistanceComparator {
    public:

        ValueEuclideanDistanceComparator(double distThreshold = 0.01, double valueThreshold = 1e-4);

        bool operator()(const Reference &lhs, const Reference &rhs) const;

    private:
        double distThreshold_;
        double valueThreshold_;
    };



}

#endif //AMOLQCPP_COMPARATORS_H
