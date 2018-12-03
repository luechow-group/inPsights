//
// Created by heuer on 03.12.18.
//

#ifndef INPSIGHTS_SPINCONNECTIONS3D_H
#define INPSIGHTS_SPINCONNECTIONS3D_H

#include <IConnection.h>
#include <ParticlesVector3D.h>

class SpinConnections3D : public IConnection {
public:
    explicit SpinConnections3D(ElectronsVector3D *electronsVector3D,
            double maxDistance = 1.6,
            double identicalThreshold = 0.01);

private:
    enum class PairType { atSamePosition, closeBy };

    double maxDistance_, identicalThreshold_; //TODO add to general settings
    std::vector<bool> atSamePositionQList_;
    std::map<std::pair<long, long>, PairType> pairIndicesMap_;

    void classifyElectrons(ElectronsVector3D *electronsVector3D);

    void createConnections(ElectronsVector3D *electronsVector3D);
};

#endif //INPSIGHTS_SPINCONNECTIONS3D_H
