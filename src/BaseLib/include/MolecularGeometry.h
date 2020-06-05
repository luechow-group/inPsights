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
