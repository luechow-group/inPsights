//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTION_H
#define AMOLQCGUI_ELECTRONCOLLECTION_H

#include "ParticleCollection.h"
#include "Electron.h"

class ElectronCollection : ParticleCollection{
    ElectronCollection(const VectorXd& positions)
            : ParticleCollection(positions)
    {
    }
};

#endif //AMOLQCGUI_ELECTRONCOLLECTION_H
