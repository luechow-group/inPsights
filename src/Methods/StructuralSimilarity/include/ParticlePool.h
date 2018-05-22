//
// Created by Michael Heuer on 08.05.18.
//

#ifndef AMOLQCPP_PARTICLEPOOL_H
#define AMOLQCPP_PARTICLEPOOL_H

#include <vector>
#include <ParticlesVector.h>
#include <MolecularGeometry.h>

using  AtomKit = std::vector<std::pair<Elements::ElementType,unsigned>>;
using  ElectronKit = std::pair<unsigned,unsigned>; // alpha, beta

namespace ParticlePool{
    void create(const AtomKit& atomKit, int charge = 0, unsigned multiplicity = 1);

    void create(const AtomsVector& atoms, int charge = 0, unsigned multiplicity = 1);

    void create(const AtomsVector& atoms, const ElectronsVector& electrons);

    AtomKit createAtomKitFromAtomsVector(const AtomsVector& atoms);

    ElectronKit createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity);

    bool isSubsetQ(const AtomsVector& atomsVector);

    unsigned numberOfAtoms();

    unsigned numberOfElectrons();

    extern AtomKit atomKit;
    extern ElectronKit electronKit;
};

#endif //AMOLQCPP_PARTICLEPOOL_H
