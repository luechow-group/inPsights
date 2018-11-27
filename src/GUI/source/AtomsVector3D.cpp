//
// Created by heuer on 24.05.17.
//

#include <AtomsVector3D.h>
#include <Atom3D.h>
#include <Bond3D.h>
#include <GuiHelper.h>

AtomsVector3D::AtomsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector) {

    /*TODO Refactor together with ElectronsVector3D*/
    std::vector<Atom3D> atoms3D;

    // Draw atoms
    for (long i = 0; i < atomsVector.numberOfEntities(); ++i) {
        atoms3D.emplace_back(Atom3D(
                root,
                GuiHelper::toQVector3D(atomsVector[i].position()),
                atomsVector.typesVector()[i]));
    }
}
