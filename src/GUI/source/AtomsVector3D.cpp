//
// Created by heuer on 24.05.17.
//

#include <AtomsVector3D.h>
#include <Particle3D.h>

AtomsVector3D::AtomsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector)
        : QEntity(root) {

    for (long i = 0; i < atomsVector.numberOfEntities(); ++i)
        atoms3D_.emplace_back(new Atom3D(root,atomsVector[i]));
}
