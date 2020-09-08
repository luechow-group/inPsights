// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MOLECULARGEOMETRY_H
#define INPSIGHTS_MOLECULARGEOMETRY_H

#include "ParticlesVector.h"
#include <list>

class MolecularGeometry{
public:
    MolecularGeometry();
    MolecularGeometry(AtomsVector atoms, ElectronsVector electrons);

    const AtomsVector& atoms() const;
    AtomsVector& atoms();

    const ElectronsVector& electrons() const;
    ElectronsVector & electrons();

    PositionsVector positions() const;

    Particle<int> operator[](long i) const;

    std::pair<bool,long> findIndexByEnumeratedType(const EnumeratedType<int> &enumeratedType) const;

    EnumeratedType<int> findEnumeratedTypeByIndex(unsigned idx) const;

    long numberOfEntities() const;

    //TODO refactor Motifs::classifyMotifs since it is the only user of the method
    std::tuple<bool,Eigen::Index> electronAtNucleusQ(long i, double threshold = 0.01) const ;

    std::list<long> coreElectronsIndices(long k, double threshold = 0.01) const;

    std::list<long> coreElectronsIndices(double threshold = 0.01) const;

    std::list<long> nonCoreElectronsIndices(double threshold = 0.01) const;

    bool operator==(const MolecularGeometry &other) const;

    bool operator!=(const MolecularGeometry &other) const;

    friend std::ostream& operator<<(std::ostream &os, const MolecularGeometry &mol) {
        os << mol.atoms() << std::endl;
        os << mol.electrons() << std::endl;
        return os;
    }

    struct Permutation {
        Eigen::PermutationMatrix<Eigen::Dynamic> nuclearPermutation, electronicPermutation;
    };

    Permutation splitAllParticlePermutation(const Eigen::PermutationMatrix<Eigen::Dynamic>& allParticlePermutation) const;

private:
    AtomsVector atoms_;
    ElectronsVector electrons_;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<>
    struct convert<MolecularGeometry> {
        static Node encode(const MolecularGeometry& rhs);
        static bool decode(const Node& node, MolecularGeometry& rhs);
    };
    Emitter& operator<< (Emitter& out, const MolecularGeometry& p);
}


#endif //INPSIGHTS_MOLECULARGEOMETRY_H
