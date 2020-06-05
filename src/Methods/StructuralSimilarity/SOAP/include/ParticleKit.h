/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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

        namespace internal {
            void createAtomKitFromAtomsVector(const AtomsVector &atoms);

            void createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity);

            void createElectronKitFromElectronsVector(const ElectronsVector &electronsVector);
        }

        ElementTypesVector toElementTypesVector();

        SpinTypesVector toSpinTypesVector();

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const AtomsVector &atomsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const AtomsVector &atomsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const ElectronsVector &electronsVector);

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const ElectronsVector &electronsVector);

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
