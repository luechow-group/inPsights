//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H
#define AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H

#include <AtomCollection.h>
#include "Abstract3dObject.h"

class MolecularGeometry3D{
public:
    MolecularGeometry3D(Qt3DCore::QEntity *root,AtomCollection atomCollection);
};

#endif //AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H
