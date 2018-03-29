//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCPP_MOLECULARGEOMETRYIMPORTER_H
#define AMOLQCPP_MOLECULARGEOMETRYIMPORTER_H

#include "Abstract3dObject.h"
#include "AtomsVector.h"

class AtomsVector3D{
public:
    AtomsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector);
};

#endif //AMOLQCPP_MOLECULARGEOMETRYIMPORTER_H
