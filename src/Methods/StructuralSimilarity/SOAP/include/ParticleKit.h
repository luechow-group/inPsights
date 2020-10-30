// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PARTICLEKIT_H
#define INPSIGHTS_PARTICLEKIT_H

#include <vector>
#include <ParticlesVector.h>
#include <MolecularGeometry.h>
#include <Eigen/Core>


using AtomKit = std::vector<std::pair<Element, unsigned>>;
using ElectronKit = std::pair<unsigned, unsigned>; // alpha, beta
using TypeKit = std::vector<std::pair<int, unsigned>>;

namespace SOAP {
    namespace ParticleKit {

        void create(const AtomKit &atomKit, const ElectronKit &electronKit);

        void createKit(const AtomKit &atomKit, const ElectronKit &electronKit);

        void create(const AtomKit &atomKit, int charge = 0, unsigned multiplicity = 1);

        void create(const AtomsVector &atoms, int charge = 0, unsigned multiplicity = 1);

        void create(const AtomsVector &atoms, const ElectronsVector &electrons);

        void create(const MolecularGeometry &molecularGeometry);

        AtomKit merge(const AtomKit& a, const AtomKit& b);

        ElectronKit merge(const ElectronKit& a, const ElectronKit& b);

        namespace internal {
            AtomKit createAtomKitFromAtomsVector(const AtomsVector &atoms);

            ElectronKit createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity);

            ElectronKit createElectronKitFromElectronsVector(const ElectronsVector &electronsVector);
        }

        ElementTypesVector toElementTypesVector();

        SpinTypesVector toSpinTypesVector();

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const AtomsVector &atomsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const AtomsVector &atomsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const ElectronsVector &electronsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const ElectronsVector &electronsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const MolecularGeometry &molecule);

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const MolecularGeometry &molecule);

        bool isSubsetQ(const AtomsVector &atomsVector);

        bool isSubsetQ(const ElectronsVector &electronsVector);

        bool isSubsetQ(const MolecularGeometry &molecularGeometry);

        bool isSameSetQ(const AtomsVector &atomsVector);

        bool isSameSetQ(const ElectronsVector &electronsVector);

        bool isSameSetQ(const MolecularGeometry &molecularGeometry);

        EnumeratedElement getEnumeratedElementByIndex(unsigned idx);

        EnumeratedSpin getEnumeratedSpinByIndex(unsigned idx);

        EnumeratedType<int> getEnumeratedTypeByIndex(unsigned idx);

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
}

#endif //INPSIGHTS_PARTICLEKIT_H
