//
// Created by heuer on 03.12.18.
//

#ifndef INPSIGHTS_SPINCLASSIFICATION_H
#define INPSIGHTS_SPINCLASSIFICATION_H

#include <ParticlesVector.h>


namespace SpinClassification{
    enum class PairType { atSamePosition, closeBy };

    using PairTypeMap = std::map<std::pair<long, long>, SpinClassification::PairType>;

    PairTypeMap classify(const ElectronsVector& electronsVector,
            double maxDistance = std::numeric_limits<double>::max(),
            double identicalThreshold = 0.01);

    bool isAtSamePositionQ(const std::map<std::pair<long, long>, PairType>& pairTypes, long i);
}


#endif //INPSIGHTS_SPINCLASSIFICATION_H
