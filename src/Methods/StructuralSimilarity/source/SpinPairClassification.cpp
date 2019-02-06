//
// Created by heuer on 03.12.18.
//

#include <SpinPairClassification.h>
#include <Metrics.h>

SpinPairClassification::PairTypeMap
SpinPairClassification::classify(const ElectronsVector& electronsVector, double maxDistance, double identicalThreshold) {

    PairTypeMap pairTypes;

    for (Eigen::Index i = 0; i < electronsVector.numberOfEntities(); ++i) {
        for (Eigen::Index j = i + 1; j < electronsVector.numberOfEntities(); ++j) {
            if (atSamePositionQ({i,j},electronsVector, identicalThreshold))
                pairTypes[{i, j}] = PairType::atSamePosition;

            else if (closeByQ({i,j}, electronsVector, maxDistance))
                pairTypes[{i, j}] = PairType::closeBy;
            else
                pairTypes[{i, j}] = PairType::distant;
        }
    }
    return pairTypes;
}

bool SpinPairClassification::atSamePositionQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
                                         const ElectronsVector &electronsVector, double identicalThreshold) {
    auto dist = Metrics::distance(
            electronsVector[pair.first].position(),
            electronsVector[pair.second].position());

    if(dist <= identicalThreshold){
        assert(electronsVector[pair.first].type() != electronsVector[pair.second].type()
               && "Antisymmetry violation! Two same spin electrons cannot have the same position.");
        return true;
    } else {
        return false;
    }
}

bool SpinPairClassification::closeByQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
              const ElectronsVector& electronsVector, double maxDistance) {
    auto dist = Metrics::distance(
            electronsVector[pair.first].position(),
            electronsVector[pair.second].position());

    return  dist <= maxDistance;
}


bool SpinPairClassification::isAtSamePositionQ(const PairTypeMap &pairTypes, Eigen::Index i) {

    for(const auto& pair : pairTypes)
        if( (pair.first.first == i || pair.first.second == i) && pair.second == PairType::atSamePosition)
            return true;
    return false;
}
