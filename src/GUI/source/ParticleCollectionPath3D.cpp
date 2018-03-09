//
// Created by Michael Heuer on 12.11.17.
//

#include <Polyline.h>
#include "ParticleCollectionPath3D.h"
#include "Electron3D.h"


ParticleCollectionPath3D::ParticleCollectionPath3D(Qt3DCore::QEntity *root,
                                                   const ElectronCollections &electronCollections,
                                                   float radius) {


    auto numberOfParticles = electronCollections[0].numberOfEntities();
    std::vector<std::vector<QVector3D>> pointsList(numberOfParticles);
    for (unsigned i = 0; i < numberOfParticles; ++i) { // iterate over particles

        for (int j = 0; j < electronCollections.numberOfEntities(); ++j) {
            auto tmp = electronCollections[j][i].position();
            pointsList[i].emplace_back(QVector3D(float(tmp(0)),float(tmp(1)),float(tmp(2))));
        }

        auto spinType = electronCollections.spinTypeCollection()[i];

        if (spinType == Spin::SpinType::alpha) {
            new Polyline(root,Spin::QColorFromSpinType(Spin::SpinType::alpha) , pointsList[i], radius);
        }
        else {
            new Polyline(root,Spin::QColorFromSpinType(Spin::SpinType::beta) , pointsList[i], radius);
        }

    }

}
