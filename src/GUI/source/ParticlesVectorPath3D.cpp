//
// Created by Michael Heuer on 12.11.17.
//

#include <Polyline.h>
#include "ParticlesVectorPath3D.h"
#include "Electron3D.h"


ParticlesVectorPath3D::ParticlesVectorPath3D(Qt3DCore::QEntity *root,
                                                   const ElectronsVectorCollection &electronsVectorCollection,
                                                   float radius) {


    auto numberOfParticles = electronsVectorCollection[0].numberOfEntities();
    std::vector<std::vector<QVector3D>> pointsList(numberOfParticles);
    for (unsigned i = 0; i < numberOfParticles; ++i) { // iterate over particles

        for (int j = 0; j < electronsVectorCollection.numberOfEntities(); ++j) {
            auto tmp = electronsVectorCollection[j][i].position();
            pointsList[i].emplace_back(QVector3D(float(tmp(0)),float(tmp(1)),float(tmp(2))));
        }

        auto spinType = electronsVectorCollection.spinTypesVector()[i];

        if (spinType == Spin::SpinType::alpha) {
            new Polyline(root,Spin::QColorFromSpinType(Spin::SpinType::alpha) , pointsList[i], radius);
        }
        else {
            new Polyline(root,Spin::QColorFromSpinType(Spin::SpinType::beta) , pointsList[i], radius);
        }

    }

}
