//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
#define AMOLQCGUI_ELECTRONCOLLECTIONPATH_H

#include "ParticleCollections.h"
#include "ElectronCollection.h"
#include "SpinTypeCollection.h"

class ElectronCollections : public ParticleCollections{
public:
    ElectronCollections();
    ElectronCollections(const std::vector<ElectronCollection>& electronCollections);
    ElectronCollections(const std::vector<ParticleCollection>& particleCollections_,
                        const SpinTypeCollection& spinTypeCollection);

private:
    SpinTypeCollection spinTypeCollection_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
