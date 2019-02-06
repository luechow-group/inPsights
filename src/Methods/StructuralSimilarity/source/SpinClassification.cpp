//
// Created by heuer on 03.12.18.
//

#include <SpinClassification.h>
#include <Metrics.h>

SpinClassification::PairTypeMap
SpinClassification::classify(const ElectronsVector& electronsVector, double maxDistance, double identicalThreshold) {

    PairTypeMap pairTypes;

    for (long i = 0; i < electronsVector.numberOfEntities(); ++i) {
        for (long j = i + 1; j < electronsVector.numberOfEntities(); ++j) {
            auto dist = Metrics::distance(
                    electronsVector[i].position(),
                    electronsVector[j].position());

            if (dist <= identicalThreshold) {
                pairTypes[{i, j}] = PairType::atSamePosition;

                assert(electronsVector[i].type() != electronsVector[j].type()
                       && "Antisymmetry violation! Types of paired electrons cannot be identical.");

            } else if (dist < maxDistance) {
                pairTypes[{i, j}] = PairType::closeBy;
            }
        }
    }
    return pairTypes;
}

bool SpinClassification::isAtSamePositionQ(
        const PairTypeMap &pairTypes, long i) {

    for(const auto& pair : pairTypes)
        if( (pair.first.first == i || pair.first.second == i) && pair.second == PairType::atSamePosition)
            return true;
    return false;
}
