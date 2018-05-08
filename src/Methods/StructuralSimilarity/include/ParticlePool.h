//
// Created by Michael Heuer on 08.05.18.
//

#ifndef AMOLQCPP_PARTICLEPOOL_H
#define AMOLQCPP_PARTICLEPOOL_H

#include <vector>
#include <ParticlesVector.h>

//TODO use std::map instead

using  AtomKit = std::vector<std::pair<Elements::ElementType,unsigned>>;
using  ElectronKit = std::pair<unsigned,unsigned>; // alpha, beta

class ParticlePool{
public:
    explicit ParticlePool(const AtomKit& atomKit, int charge = 0, unsigned multiplicity = 1)
            : atomKit_(atomKit),
              electronKit_(createElectronKitFromAtomKit(atomKit, charge, multiplicity))
    {}

    explicit ParticlePool(const AtomsVector& atoms, int charge = 0, unsigned multiplicity = 1)
            : ParticlePool(createAtomKitFromAtomsVector(atoms), charge, multiplicity)
    {}

    ParticlePool(const AtomsVector& atoms, const ElectronsVector& electrons)
            : atomKit_(createAtomKitFromAtomsVector(atoms)),
              electronKit_{electrons.countTypeOccurence(Spins::SpinType::alpha),
                           electrons.countTypeOccurence(Spins::SpinType::beta)}
    {}

    AtomKit createAtomKitFromAtomsVector(const AtomsVector& atoms){

        AtomKit atomKit;
        std::vector<Elements::ElementType> elementsPresent;

        for (int i = 0; i < atoms.numberOfEntities(); ++i) {
            auto element = atoms[i].type();
            if (std::find(elementsPresent.begin(), elementsPresent.end(), element) ==elementsPresent.end())
                elementsPresent.emplace_back(element);
        }

        std::sort(elementsPresent.begin(), elementsPresent.end());

        for (const auto& element : elementsPresent){
            atomKit.emplace_back<std::pair<Elements::ElementType,unsigned>>({element,atoms.countTypeOccurence(element)});
        }

        return atomKit;
    };

    ElectronKit createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity){
        unsigned numberOfElectrons = 0;

        for (auto const &elemenTypeNumberPair : atomKit_) {
            auto elementType = elemenTypeNumberPair.first;
            auto numberOfAtoms = elemenTypeNumberPair.second;

            numberOfElectrons += Elements::ElementInfo::Z(elementType)*numberOfAtoms;
        }
        numberOfElectrons -= charge;

        unsigned numberOfUnpairedElectrons = multiplicity-1;

        assert((numberOfElectrons - numberOfUnpairedElectrons)%2 == 0
               && "The number of electron pairs must be a multiple of 2.");
        unsigned numberOfElectronPairs = (numberOfElectrons - numberOfUnpairedElectrons)/2;

        // use amolqc convention
        unsigned numberOfAlphaElectrons = numberOfUnpairedElectrons + numberOfElectronPairs;
        unsigned numberOfBetaElectrons = numberOfElectronPairs;


        return {numberOfAlphaElectrons,numberOfBetaElectrons};
    };

    bool isSubsetQ(const AtomsVector& atomsVector){
        bool isSubsetQ = true;

        for (auto const &elemenTypeNumberPair : atomKit_) {
            auto elementType = elemenTypeNumberPair.first;
            auto maxNumberOfAtoms = elemenTypeNumberPair.second;

            unsigned short count = 0;
            for (int j = 0; j < atomsVector.numberOfEntities(); ++j) {

                if(atomsVector[j].type() == elementType)
                    count++;
                if( count > maxNumberOfAtoms){
                    isSubsetQ = false;
                    goto outsideNestedLoop;
                }
            }
        }

        outsideNestedLoop:
        return isSubsetQ;
    }

    const AtomKit& atomKit() const {
        return atomKit_;
    }

    const ElectronKit& electronKit() const {
        return electronKit_;
    }


private:
    AtomKit atomKit_;
    ElectronKit electronKit_;
};

#endif //AMOLQCPP_PARTICLEPOOL_H
