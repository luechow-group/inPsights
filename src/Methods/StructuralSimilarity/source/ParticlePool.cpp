//
// Created by Michael Heuer on 15.05.18.
//

#include "ParticlePool.h"

AtomKit ParticlePool::atomKit_ = {};

ElectronKit ParticlePool::electronKit_ = {0,0};

ParticlePool::ParticlePool(const AtomKit& atomKit, int charge, unsigned multiplicity) {
    atomKit_ = atomKit;

    //TODO refactor
    auto ekit = createElectronKitFromAtomKit(atomKit, charge, multiplicity);
    electronKit_.first = ekit.first;
    electronKit_.second = ekit.second;
}

ParticlePool::ParticlePool(const AtomsVector& atoms, int charge, unsigned multiplicity)
        : ParticlePool(createAtomKitFromAtomsVector(atoms), charge, multiplicity)
{}

ParticlePool::ParticlePool(const AtomsVector& atoms, const ElectronsVector& electrons) {
    atomKit_ = createAtomKitFromAtomsVector(atoms);

    electronKit_.first = electrons.countTypeOccurence(Spins::SpinType::alpha);
    electronKit_.second = electrons.countTypeOccurence(Spins::SpinType::beta);
}

AtomKit ParticlePool::createAtomKitFromAtomsVector(const AtomsVector& atoms) const {

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

ElectronKit ParticlePool::createElectronKitFromAtomKit(const AtomKit &atomKit,
                                                       int charge, unsigned multiplicity) const {
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

bool ParticlePool::isSubsetQ(const AtomsVector& atomsVector) const {
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

const AtomKit & ParticlePool::atomKit() {
    return atomKit_;
}

const ElectronKit& ParticlePool::electronKit() {
    return electronKit_;
}
