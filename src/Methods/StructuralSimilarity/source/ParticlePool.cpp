//
// Created by Michael Heuer on 15.05.18.
//

#include "ParticlePool.h"

namespace ParticleKit {
    AtomKit atomKit = {};

    ElectronKit electronKit = {0, 0};

    void create(const AtomKit &atomKit, const ElectronKit &electronKit) {
        ParticleKit::atomKit = atomKit;
        ParticleKit::electronKit.first = electronKit.first;
        ParticleKit::electronKit.second = electronKit.second;
    }

    void create(const AtomKit &atomKit, int charge, unsigned multiplicity) {
        ParticleKit::atomKit = atomKit;

        createElectronKitFromAtomKit(atomKit, charge, multiplicity);
    }

    void create(const AtomsVector &atoms, int charge, unsigned multiplicity){
        createAtomKitFromAtomsVector(atoms);
        create(ParticleKit::atomKit, charge, multiplicity);
    }

    void create(const AtomsVector &atoms, const ElectronsVector &electrons) {
        createAtomKitFromAtomsVector(atoms);
        createElectronKitFromElectronsVector(electrons);
    }

    namespace {
        void createAtomKitFromAtomsVector(const AtomsVector &atoms) {

            AtomKit newAtomKit;
            std::vector<Elements::ElementType> elementsPresent;

            for (int i = 0; i < atoms.numberOfEntities(); ++i) {
                auto element = atoms[i].type();
                if (std::find(elementsPresent.begin(), elementsPresent.end(), element) == elementsPresent.end())
                    elementsPresent.emplace_back(element);
            }

            std::sort(elementsPresent.begin(), elementsPresent.end());

            for (const auto &element : elementsPresent) {
                newAtomKit.emplace_back<std::pair<Elements::ElementType, unsigned>>(
                        {element, atoms.countTypeOccurence(element)});
            }

            ParticleKit::atomKit = newAtomKit;
        };


        void createElectronKitFromElectronsVector(const ElectronsVector &electronsVector) {

            electronKit.first = electronsVector.countTypeOccurence(Spins::SpinType::alpha);
            electronKit.second = electronsVector.countTypeOccurence(Spins::SpinType::beta);
        }

        void createElectronKitFromAtomKit(const AtomKit &atomKit,
                                          int charge, unsigned multiplicity) {
            unsigned numberOfElectrons = 0;

            for (auto const &elemenTypeNumberPair : atomKit) {
                auto elementType = elemenTypeNumberPair.first;
                auto numberOfAtoms = elemenTypeNumberPair.second;

                numberOfElectrons += Elements::ElementInfo::Z(elementType) * numberOfAtoms;
            }
            numberOfElectrons -= charge;

            unsigned numberOfUnpairedElectrons = multiplicity - 1;

            assert((numberOfElectrons - numberOfUnpairedElectrons) % 2 == 0
                   && "The number of electron pairs must be a multiple of 2.");
            unsigned numberOfElectronPairs = (numberOfElectrons - numberOfUnpairedElectrons) / 2;

            // use amolqc convention
            unsigned numberOfAlphaElectrons = numberOfUnpairedElectrons + numberOfElectronPairs;
            unsigned numberOfBetaElectrons = numberOfElectronPairs;

            ParticleKit::electronKit.first = numberOfAlphaElectrons;
            ParticleKit::electronKit.second = numberOfBetaElectrons;
        };
    }

    bool isSubsetQ(const AtomsVector &atomsVector) {
        bool isSubsetQ = true;

        for (auto const &elemenTypeNumberPair : atomKit) {
            auto elementType = elemenTypeNumberPair.first;
            auto maxNumberOfAtoms = elemenTypeNumberPair.second;

            unsigned short count = 0;
            for (int j = 0; j < atomsVector.numberOfEntities(); ++j) {

                if (atomsVector[j].type() == elementType)
                    count++;
                if (count > maxNumberOfAtoms) {
                    isSubsetQ = false;
                    goto outsideNestedLoop;
                }
            }
        }

        outsideNestedLoop:
        return isSubsetQ;
    }

    unsigned numberOfAtoms() {
        unsigned sum = 0;

        for (auto const &elementTypeNumberPair :  atomKit)
            sum += elementTypeNumberPair.second;

        return sum;
    }

    unsigned numberOfElectrons() {
        return electronKit.first + electronKit.second;
    }
}
