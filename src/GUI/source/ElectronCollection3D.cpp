//
// Created by Michael Heuer on 12.11.17.
//

#include "ElectronCollection3D.h"
#include "Electron3D.h"
#include "Bond3D.h"

ElectronCollection3D::ElectronCollection3D(Qt3DCore::QEntity *root, const ElectronCollection &electronCollection) {

    std::vector<Electron3D> electrons3D;

    // Draw electrons
    for (long i = 0; i < electronCollection.numberOfParticles(); ++i) {
        Eigen::Vector3d vec= electronCollection[i].position();
        electrons3D.emplace_back(Electron3D(root, QVector3D(float(vec[0]),float(vec[1]),float(vec[2])),
                                            electronCollection.spinType(i)));
    }

}
