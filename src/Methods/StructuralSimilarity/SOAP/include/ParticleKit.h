//
// Created by Michael Heuer on 08.05.18.
//

#ifndef INPSIGHTS_PARTICLEKIT_H
#define INPSIGHTS_PARTICLEKIT_H

#include <vector>
#include <ParticlesVector.h>
#include <MolecularGeometry.h>

using AtomKit = std::vector<std::pair<Element,unsigned>>;
using ElectronKit = std::pair<unsigned,unsigned>; // alpha, beta
using TypeKit = std::vector<std::pair<int,unsigned>>;

namespace ParticleKit{

    void create(const AtomKit& atomKit, const ElectronKit& electronKit);

    void createKit(const AtomKit &atomKit, const ElectronKit &electronKit);

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

    NumberedElement getNumberedElementByIndex(unsigned idx);

    NumberedSpin getNumberedSpinByIndex(unsigned idx);

    NumberedType<int> getNumberedTypeByIndex(unsigned idx);

    unsigned numberOfElementTypes();

    unsigned numberOfSpinTypes();

    unsigned numberOfTypes();

    unsigned numberOfAtoms();

    unsigned numberOfElectrons();

    unsigned numberOfParticles();

    std::string toString();

    extern AtomKit atomKit;
    extern ElectronKit electronKit;
    extern TypeKit kit;
};

#endif //INPSIGHTS_PARTICLEKIT_H
