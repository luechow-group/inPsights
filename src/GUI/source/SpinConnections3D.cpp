//
// Created by heuer on 03.12.18.
//

#include <SpinConnections3D.h>
#include <Metrics.h>

SpinConnections3D::SpinConnections3D(ElectronsVector3D *electronsVector3D,
                                     double identicalThreshold,
                                     double maxDistance)
        : IConnection(electronsVector3D->connections_),
          identicalThreshold_(identicalThreshold),
          maxDistance_(maxDistance),
          atSamePositionQList_(size_t(electronsVector3D->numberOfEntities()), false),
          pairIndicesMap_() {

    classifyElectrons(electronsVector3D);
}

void SpinConnections3D::classifyElectrons(ElectronsVector3D *electronsVector3D) {
    for (long i = 0; i < electronsVector3D->numberOfEntities(); ++i) {
        for (long j = i + 1; j < electronsVector3D->numberOfEntities(); ++j) {
            auto dist = Metrics::distance(
                    electronsVector3D->operator[](i).position(),
                    electronsVector3D->operator[](j).position());

            if (dist <= identicalThreshold_) {
                atSamePositionQList_[i] = true;
                atSamePositionQList_[j] = true;
                pairIndicesMap_[{i, j}] = PairType::atSamePosition;

                assert(electronsVector3D->operator[](i).type() != electronsVector3D->operator[](j).type()
                       && "Antisymmetry violation! Types of paired electrons cannot be identical.");

            } else if (dist < maxDistance_) {
                pairIndicesMap_[{i, j}] = PairType::closeBy;
            }
        }
    }
}

