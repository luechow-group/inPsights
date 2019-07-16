//
// Created by Michael Heuer on 2019-04-23.
//

#include "IClusterer.h"

IClusterer::Type IClusterer::typeFromString(const std::string& clustererName) {
    if(clustererName == "BestMatchDistanceIdentityClusterer")
        return Type::BestMatchDistanceIdentityClusterer;
    else if(clustererName == "BestMatchDistanceSimilarityClusterer")
        return Type::BestMatchDistanceSimilarityClusterer;
    else if(clustererName == "BestMatchDistanceDensityBasedClusterer")
        return Type::BestMatchDistanceDensityBasedClusterer;
    else if(clustererName == "ReferencePositionsClusterer")
        return Type::ReferencePositionsClusterer;
    else if(clustererName == "BestMatchSOAPSimilarityClusterer")
        return Type::BestMatchSOAPSimilarityClusterer;
    else
        return Type::invalid;
};
