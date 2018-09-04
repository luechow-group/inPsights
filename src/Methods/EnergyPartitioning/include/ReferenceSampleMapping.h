//
// Created by Michael Heuer on 31.08.18.
//

#ifndef AMOLQCPP_REFERENCESAMPLEMAPPING_H
#define AMOLQCPP_REFERENCESAMPLEMAPPING_H

#include <map>
#include "SampleData.h"
#include "Comparators.h"

using RefSamplePair = std::pair<Reference,Sample>;//TODO used?

class ReferenceSampleMapping {
public:
    std::multimap<Reference,Sample,Comparators::ValueEuclideanDistanceComparator> map;
};

#endif //AMOLQCPP_REFERENCESAMPLEMAPPING_H
