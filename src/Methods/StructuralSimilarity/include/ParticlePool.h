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

        //TODO MAKE MONOSTATE
class ParticlePool{
public:
    explicit ParticlePool(const AtomKit& atomKit, int charge = 0, unsigned multiplicity = 1);

    explicit ParticlePool(const AtomsVector& atoms, int charge = 0, unsigned multiplicity = 1);

    ParticlePool(const AtomsVector& atoms, const ElectronsVector& electrons);

    AtomKit createAtomKitFromAtomsVector(const AtomsVector& atoms) const;

    ElectronKit createElectronKitFromAtomKit(const AtomKit &atomKit, int charge, unsigned multiplicity) const;

    bool isSubsetQ(const AtomsVector& atomsVector) const;

    static const AtomKit& atomKit();

    static const ElectronKit& electronKit();

private:
    static AtomKit atomKit_;
    static ElectronKit electronKit_;
};

#endif //AMOLQCPP_PARTICLEPOOL_H
