//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H
#define AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H

#include "Abstract3dObject.h"
#include "AtomsVector.h"

class AtomsVector3D{
public:
    AtomsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector);
};

#endif //AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H
