//
// Created by Michael Heuer on 2019-04-12.
//

#ifndef INPSIGHTS_ICLUSTERER_H
#define INPSIGHTS_ICLUSTERER_H

#include <Group.h>

class IClusterer{
public:
    enum class Type {
        IdentityClusterer,
        DistanceClusterer,
        BestMatchDistanceDensityBasedClusterer,
        BestMatchSOAPSimilarityClusterer,
        ReferencePositionsClusterer,
        invalid
    };

    virtual void cluster(Group& group) = 0;

    static Type typeFromString(const std::string& clustererName);
};

#endif //INPSIGHTS_ICLUSTERER_H
