//
// Created by Michael Heuer on 08.05.18.
//

#ifndef AMOLQCPP_PARTICLEKIT_H
#define AMOLQCPP_PARTICLEKIT_H

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

    void create(const MolecularGeometry& molecularGeometry);

    namespace {
        void createAtomKitFromAtomsVector(const AtomsVector& atoms);

        void createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity);

        void createElectronKitFromElectronsVector(const ElectronsVector &electronsVector);
    }

    bool isSubsetQ(const AtomsVector& atomsVector);

    bool isSubsetQ(const ElectronsVector& electronsVector);

    bool isSubsetQ(const MolecularGeometry& molecularGeometry);

    unsigned numberOfTypes();

    unsigned numberOfAtoms();

    unsigned numberOfElectrons();

    unsigned numberOfParticles();

    extern AtomKit atomKit;
    extern ElectronKit electronKit;
};

#endif //AMOLQCPP_PARTICLEKIT_H
