//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H
#define AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H

#include "Abstract3dObject.h"
#include "AtomCollection.h"

class AtomCollection3D{
public:
    AtomCollection3D(Qt3DCore::QEntity *root, const AtomCollection &atomCollection);
};

#endif //AMOLQCGUI_MOLECULARGEOMETRYIMPORTER_H
