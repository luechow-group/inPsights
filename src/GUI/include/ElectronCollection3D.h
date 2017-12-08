//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTION3D_H
#define AMOLQCGUI_ELECTRONCOLLECTION3D_H

#include "Abstract3dObject.h"
#include "ElectronCollection.h"

class ElectronCollection3D{
public:
    ElectronCollection3D(Qt3DCore::QEntity *root, const ElectronCollection &electonCollection,
                             bool showIndicesQ = false);
};

#endif //AMOLQCGUI_ELECTRONCOLLECTION3D_H
