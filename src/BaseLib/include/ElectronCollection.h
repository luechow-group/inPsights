//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTION_H
#define AMOLQCGUI_ELECTRONCOLLECTION_H

#include "ParticleCollection.h"
#include "SpinTypeCollection.h"
#include "Electron.h"

class ElectronCollection : public ParticleCollection, public SpinTypeCollection{
public:
    explicit ElectronCollection(const VectorXd& positions);
    explicit ElectronCollection(const VectorXd& positions, const VectorXi& spinTypes);

    Electron electron(long i);
    Spin::SpinType spinType(long i);

private:
    //SpinTypeCollection spinTypes_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTION_H
