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

namespace ParticleKit{

    void create(const AtomKit& atomKit, const ElectronKit& electronKit);

    void create(const AtomKit& atomKit, int charge = 0, unsigned multiplicity = 1);

    void create(const AtomsVector& atoms, int charge = 0, unsigned multiplicity = 1);

    void create(const AtomsVector& atoms, const ElectronsVector& electrons);

    namespace {
        void createAtomKitFromAtomsVector(const AtomsVector& atoms);

        void createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity);

        void createElectronKitFromElectronsVector(const ElectronsVector &electronsVector);
    }

    bool isSubsetQ(const AtomsVector& atomsVector);

    unsigned numberOfAtoms();

    unsigned numberOfElectrons();

    extern AtomKit atomKit;
    extern ElectronKit electronKit;
};

#endif //AMOLQCPP_PARTICLEPOOL_H
