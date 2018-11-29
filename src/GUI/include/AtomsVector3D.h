//
// Created by heuer on 23.05.17.
//

#ifndef INPSIGHTS_MOLECULARGEOMETRYIMPORTER_H
#define INPSIGHTS_MOLECULARGEOMETRYIMPORTER_H

#include <Qt3DCore/QEntity>
#include <ParticlesVector.h>
#include <Particle3D.h>

class AtomsVector3D : public Qt3DCore::QEntity{
public:
    AtomsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector);

    std::vector<Atom3D*> atoms3D_;
};

#endif //INPSIGHTS_MOLECULARGEOMETRYIMPORTER_H
