// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <ParticleKit.h>
#include "ElementInfo.h"

namespace SOAP {
    namespace ParticleKit {
        AtomKit atomKit = {};

        ElectronKit electronKit = {0, 0};
        TypeKit kit = {};

        void create(const AtomKit &newAtomKit, const ElectronKit &newElectronKit) {
            SOAP::ParticleKit::atomKit = newAtomKit;
            SOAP::ParticleKit::electronKit.first = newElectronKit.first;
            SOAP::ParticleKit::electronKit.second = newElectronKit.second;

            SOAP::ParticleKit::createKit(newAtomKit, electronKit);
        }

        void createKit(const AtomKit &newAtomKit, const ElectronKit &newElectronKit) {
            SOAP::ParticleKit::kit = {};

            for (const auto &[type, numberOfAtoms] : newAtomKit)
                SOAP::ParticleKit::kit.push_back({Elements::elementToInt(type), numberOfAtoms});

            auto[numberOfAlphaSpins, numberOfBetaSpins] = newElectronKit;

            if (numberOfAlphaSpins > 0)
                SOAP::ParticleKit::kit.push_back({Spins::spinToInt(Spin::alpha), numberOfAlphaSpins});

            if (numberOfBetaSpins > 0)
                SOAP::ParticleKit::kit.push_back({Spins::spinToInt(Spin::beta), numberOfBetaSpins});
        }

        void create(const AtomKit &newAtomKit, int charge, unsigned multiplicity) {
            auto newElectronKit =
                    SOAP::ParticleKit::internal::createElectronKitFromAtomKit(newAtomKit, charge, multiplicity);
            create(newAtomKit, newElectronKit);
        }

        void create(const AtomsVector &atoms, int charge, unsigned multiplicity) {
            auto newAtomKit = SOAP::ParticleKit::internal::createAtomKitFromAtomsVector(atoms);
            create(newAtomKit, charge, multiplicity);
        }

        void create(const AtomsVector &atoms, const ElectronsVector &electrons) {
            auto newAtomKit = SOAP::ParticleKit::internal::createAtomKitFromAtomsVector(atoms);
            auto newElectronKit = SOAP::ParticleKit::internal::createElectronKitFromElectronsVector(electrons);
            create(newAtomKit, newElectronKit);
        }

        void create(const MolecularGeometry &molecularGeometry) {
            create(molecularGeometry.atoms(), molecularGeometry.electrons());
        }

        ElectronKit merge(const ElectronKit& a, const ElectronKit& b) {
            auto mergedKit = a;

            if(b.first > mergedKit.first)
                mergedKit.first = b.first;
            if(b.second > mergedKit.second)
                mergedKit.second = b.second;

            return mergedKit;
        }

        AtomKit merge(const AtomKit& a, const AtomKit& b){
            auto mergedKit = a;

            // Check for all members of the second kit, if they are already of the first kit
            // and check if occurence number must be increased
            for(auto bKitMember : b) {
                auto it = std::find_if(mergedKit.begin(), mergedKit.end(),
                                        [&bKitMember](const std::pair<Element, unsigned>& aKitMember) {
                                            return aKitMember.first == bKitMember.first;
                                        });
                // if already present
                if(it != mergedKit.end()) {
                    // and number in kit B is larger
                    if (it->second < bKitMember.second)
                        it->second = bKitMember.second;
                }
                // if not present
                else
                    mergedKit.push_back(bKitMember);
            }
            return mergedKit;
        }

        namespace internal {
            AtomKit createAtomKitFromAtomsVector(const AtomsVector &atoms) {

                AtomKit newAtomKit;
                std::vector<Element> elementsPresent;

                for (int i = 0; i < atoms.numberOfEntities(); ++i) {
                    auto element = atoms[i].type();
                    if (std::find(elementsPresent.begin(), elementsPresent.end(), element) == elementsPresent.end())
                        elementsPresent.emplace_back(element);
                }

                std::sort(elementsPresent.begin(), elementsPresent.end());

                for (const auto &element : elementsPresent) {
                    newAtomKit.emplace_back<std::pair<Element, unsigned>>(
                            {element, atoms.typesVector().countOccurence(element)});
                }

                return newAtomKit;
            }

            ElectronKit createElectronKitFromElectronsVector(const ElectronsVector &electronsVector) {
                return {electronsVector.typesVector().countOccurence(Spin::alpha),
                        electronsVector.typesVector().countOccurence(Spin::beta)};
            }

            ElectronKit createElectronKitFromAtomKit(const AtomKit &atomKit,
                                              int charge, unsigned multiplicity) {
                unsigned numberOfElectrons = 0;

                for (auto const &[elementType, numberOfAtoms] : atomKit)
                    numberOfElectrons += Elements::ElementInfo::Z(elementType) * numberOfAtoms;

                numberOfElectrons -= charge;

                unsigned numberOfUnpairedElectrons = multiplicity - 1;

                // check if the number of electrons and multiplicity is correct.
                assert((numberOfElectrons - numberOfUnpairedElectrons) % 2 == 0
                       && "The number of electron pairs must be a multiple of 2.");
                unsigned numberOfElectronPairs = (numberOfElectrons - numberOfUnpairedElectrons) / 2;

                // use amolqc convention
                unsigned numberOfAlphaElectrons = numberOfUnpairedElectrons + numberOfElectronPairs;
                unsigned numberOfBetaElectrons = numberOfElectronPairs;

                return {numberOfAlphaElectrons, numberOfBetaElectrons};
            }
        }

        SpinTypesVector toSpinTypesVector() {
            return {SOAP::ParticleKit::electronKit.first, SOAP::ParticleKit::electronKit.second};
        }

        ElementTypesVector toElementTypesVector() {
            ElementTypesVector elementTypesVector;
            for (const auto &elementCount : atomKit)
                for (unsigned i = 0; i < elementCount.second; ++i)
                    elementTypesVector.append(elementCount.first);
            return elementTypesVector;
        }


        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const ElectronsVector &electronsVector) {
            assert(SOAP::ParticleKit::isSubsetQ(electronsVector));

            Eigen::VectorXi indices(electronsVector.numberOfEntities());

            int index = 0;
            for (long i = 0; i < electronsVector.numberOfEntities(); ++i)
                if (electronsVector[i].type() == Spin::alpha) {
                    indices[index] = i;
                    index++;
                }

            for (long i = 0; i < electronsVector.numberOfEntities(); ++i)
                if (electronsVector[i].type() == Spin::beta) {
                    indices[index] = i;
                    index++;
                }

            return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
        }

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const ElectronsVector &electronsVector) {
            return fromKitPermutation(electronsVector).inverse();
        }

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const AtomsVector &atomsVector) {
            assert(SOAP::ParticleKit::isSubsetQ(atomsVector));

            Eigen::VectorXi indices(atomsVector.numberOfEntities());
            auto typeCounts = atomsVector.typesVector().countTypes();

            int index = 0;
            for (const auto &[type, count] : typeCounts)
                for (long i = 0; i < atomsVector.numberOfEntities(); ++i)
                    if (atomsVector[i].type() == type) {
                        indices[index] = i;
                        index++;
                    }

            return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
        }

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const AtomsVector &atomsVector) {
            return fromKitPermutation(atomsVector).inverse();
        }

        Eigen::PermutationMatrix<Eigen::Dynamic> fromKitPermutation(const MolecularGeometry &molecule) {
            assert(SOAP::ParticleKit::isSameSetQ(molecule));
            Eigen::VectorXi indices(molecule.numberOfEntities());

            auto nuclearIndices = fromKitPermutation(molecule.atoms()).indices();
            auto M = nuclearIndices.size();

            auto electronicIndices = fromKitPermutation(molecule.electrons()).indices();
            auto N = electronicIndices.size();

            indices.head(M) = nuclearIndices;
            indices.tail(N) = electronicIndices.array() + M;

            return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
        }

        Eigen::PermutationMatrix<Eigen::Dynamic> toKitPermutation(const MolecularGeometry &molecule) {
            return fromKitPermutation(molecule).inverse();
        }

        bool isSubsetQ(const AtomsVector &atomsVector) {
            auto countedTypes = atomsVector.typesVector().countTypes();
            for (auto typeCount: countedTypes) {

                auto it = std::find_if(SOAP::ParticleKit::atomKit.begin(), SOAP::ParticleKit::atomKit.end(),
                                       [typeCount](const std::pair<Element, unsigned> &element) {
                                           return element.first == typeCount.type_;
                                       });
                if ((it == SOAP::ParticleKit::atomKit.end()) || ((*it).second < typeCount.number_))
                    return false;
            }
            return true;
        }

        bool isSubsetQ(const ElectronsVector &electronsVector) {
            return (electronsVector.typesVector().countOccurence(Spin::alpha) <= SOAP::ParticleKit::electronKit.first)
                   && (electronsVector.typesVector().countOccurence(Spin::beta) <= SOAP::ParticleKit::electronKit.second);
        }

        bool isSubsetQ(const MolecularGeometry &molecularGeometry) {
            return isSubsetQ(molecularGeometry.atoms()) && isSubsetQ(molecularGeometry.electrons());
        }

        bool isSameSetQ(const AtomsVector &atomsVector) {
            auto countedTypes = atomsVector.typesVector().countTypes();
            for (auto typeCount: countedTypes) {

                auto it = std::find_if(SOAP::ParticleKit::atomKit.begin(), SOAP::ParticleKit::atomKit.end(),
                                       [typeCount](const std::pair<Element, unsigned> &element) {
                                           return element.first == typeCount.type_;
                                       });
                if ((it == SOAP::ParticleKit::atomKit.end()) || ((*it).second != typeCount.number_))
                    return false;
            }
            return true;
        }

        bool isSameSetQ(const ElectronsVector &electronsVector) {
            return (electronsVector.typesVector().countOccurence(Spin::alpha) == SOAP::ParticleKit::electronKit.first)
            && (electronsVector.typesVector().countOccurence(Spin::beta) == SOAP::ParticleKit::electronKit.second);
        }

        bool isSameSetQ(const MolecularGeometry &molecularGeometry) {
            return isSameSetQ(molecularGeometry.atoms()) && isSameSetQ(molecularGeometry.electrons());
        }

        unsigned numberOfTypes() {
            return unsigned(kit.size());
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

        unsigned numberOfParticles() {
            return numberOfAtoms() + numberOfElectrons();
        }

        unsigned numberOfElementTypes() {
            return unsigned(atomKit.size());
        }

        unsigned numberOfSpinTypes() {
            unsigned numberOfSpinTypes = 0;
            if (electronKit.first > 0) // alpha
                numberOfSpinTypes += 1;
            if (electronKit.second > 0) // beta
                numberOfSpinTypes += 1;

            return numberOfSpinTypes;
        }

        EnumeratedType<int> getEnumeratedTypeByIndex(unsigned idx) {
            assert(idx < SOAP::ParticleKit::numberOfParticles());

            unsigned count = 0;
            for (auto &typeNumberPair : SOAP::ParticleKit::kit) {
                if (idx >= count + typeNumberPair.second) {
                    count += typeNumberPair.second;
                } else {
                    return {typeNumberPair.first, idx - count};
                }
            }
            return {0, 0}; //TODO undefined behavior
        }

        EnumeratedElement getEnumeratedElementByIndex(unsigned idx) {
            assert(idx < numberOfAtoms());

            unsigned count = 0;
            for (auto &typeNumberPair : SOAP::ParticleKit::atomKit) {
                if (idx >= count + typeNumberPair.second) {
                    count += typeNumberPair.second;
                } else {
                    return {typeNumberPair.first, idx - count};
                }
            }
            return {Element::none, 0}; //TODO undefined behavior
        }

        EnumeratedSpin getEnumeratedSpinByIndex(unsigned idx) {
            assert(idx < numberOfElectrons());

            if (idx < electronKit.first)
                return {Spin::alpha, idx};
            else
                return {Spin::beta, idx - electronKit.first};
        }

        std::string toString() {
            std::stringstream ss;

            ss << "ParticleKit:" << std::endl
               << "------------" << std::endl;

            for (auto &typeNumberPair : SOAP::ParticleKit::atomKit) {
                ss << typeNumberPair.second << "*" << typeNumberPair.first << ", ";
            }
            ss << electronKit.first << "*e" << Spin::alpha << ", "
               << electronKit.second << "*e" << Spin::beta << std::endl;
            return ss.str();
        }
    }
}