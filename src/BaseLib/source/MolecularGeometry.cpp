// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <utility>
#include "MolecularGeometry.h"
#include <Metrics.h>
#include <algorithm>
#include <numeric>

MolecularGeometry::MolecularGeometry()
        : atoms_(),electrons_()
{};

MolecularGeometry::MolecularGeometry(AtomsVector atoms, ElectronsVector electrons)
        : atoms_(std::move(atoms)),
          electrons_(std::move(electrons)) {
};

Particle<int> MolecularGeometry::operator[](long i) const {
    auto numberOfAtoms = atoms_.numberOfEntities();
    if (i < numberOfAtoms)
        return {atoms_[i].position(),
                int(atoms_[i].type())};
    else
        return {electrons_[i-numberOfAtoms].position(),
                int(electrons_[i-numberOfAtoms].type())};
}

const AtomsVector& MolecularGeometry::atoms() const { return atoms_; }

AtomsVector& MolecularGeometry::atoms() { return atoms_; }

const ElectronsVector& MolecularGeometry::electrons() const { return electrons_; }

ElectronsVector & MolecularGeometry::electrons() { return electrons_; }

PositionsVector MolecularGeometry::positions() const {
    assert(atoms().positionsVector().entityLength() == electrons().positionsVector().entityLength());

    auto entityLength = atoms().positionsVector().entityLength();

    Eigen::VectorXd positions(entityLength * numberOfEntities());

    positions.head(atoms().numberOfEntities() * entityLength) = atoms().positionsVector().asEigenVector();
    positions.tail(electrons().numberOfEntities() * entityLength) = electrons().positionsVector().asEigenVector();

    return PositionsVector(positions);
}

long MolecularGeometry::numberOfEntities() const {
    return atoms_.numberOfEntities() + electrons_.numberOfEntities();
}


EnumeratedType<int> MolecularGeometry::findEnumeratedTypeByIndex(unsigned idx) const {
    assert(idx < numberOfEntities() && "The index cannot be greater than the number of particles - 1");

    auto M = atoms().numberOfEntities();
    if(idx < M) {
        return atoms().typesVector().getEnumeratedTypeByIndex(idx).toIntType();
    } else {
        return electrons().typesVector().getEnumeratedTypeByIndex(idx-M).toIntType();
    }
}

std::pair<bool,long> MolecularGeometry::findIndexByEnumeratedType(const EnumeratedType<int> &enumeratedType) const {
    if(enumeratedType.type_ >= int(Spins::first()) && enumeratedType.type_ <= int(Spins::last())) {
        auto boolIdx = electrons().typesVector().findIndexOfEnumeratedType(
                EnumeratedSpin(Spins::spinFromInt(enumeratedType.type_), enumeratedType.number_));

        boolIdx.second += atoms().numberOfEntities(); // TODO is this the way it should be?
        return boolIdx;
    } else if(enumeratedType.type_ >= int(Elements::first()) && enumeratedType.type_ <= int(Elements::last())) {
        return atoms().typesVector().findIndexOfEnumeratedType(
                EnumeratedElement(Elements::elementFromInt(enumeratedType.type_), enumeratedType.number_));
    } else {
        return std::pair<bool,long>(false,0);
    }
}

std::tuple<bool,Eigen::Index> MolecularGeometry::electronAtNucleusQ(long i, double threshold) const {
    for (Eigen::Index k = 0; k < atoms_.numberOfEntities(); ++k)
        if(Metrics::distance(electrons_[i].position(), atoms_[k].position()) <= threshold)
            return std::pair<bool,long>(true, k);

    return std::pair<bool,long>(false,0);
}

//TODO REPLACE

std::list<long> MolecularGeometry::coreElectronsIndices(long k, double threshold) const {
    std::list<long> indices{};

    for (long i = 0; i < electrons().numberOfEntities(); ++i)
        if (Metrics::distance(electrons_[i].position(), atoms_[k].position()) <= threshold)
            indices.emplace_back(i);

    return indices;
}

std::list<long> MolecularGeometry::coreElectronsIndices(double threshold) const {
    std::list<long> indices{};

    for (long k = 0; k < atoms_.numberOfEntities(); ++k)
        indices.splice(indices.end(), coreElectronsIndices(k, threshold));

    indices.sort();
    indices.erase(std::unique( indices.begin(), indices.end() ), indices.end()); // erase duplicates
    return indices;
}

std::list<long> MolecularGeometry::nonCoreElectronsIndices(double threshold) const {
    std::list<long> diff{}, indices = std::list<long>(size_t(electrons().numberOfEntities()));
    std::iota(indices.begin(), indices.end(), 0);

    auto coreIndices = coreElectronsIndices(threshold);
    std::set_difference(indices.begin(), indices.end(), coreIndices.begin(), coreIndices.end(),
                        std::inserter(diff, diff.begin()));
    return diff;
}

bool MolecularGeometry::operator==(const MolecularGeometry &other) const {
    return (atoms() == other.atoms())
           && (electrons() == other.electrons());
}

bool MolecularGeometry::operator!=(const MolecularGeometry &other) const {
    return !(*this == other);
}

MolecularGeometry::Permutation MolecularGeometry::splitAllParticlePermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &allParticlePermutation) const {

    assert(allParticlePermutation.indices().size() == numberOfEntities());

    auto M = atoms().numberOfEntities();
    auto N = electrons().numberOfEntities();

    Eigen::VectorXi permutedNuclearIndices(M), permutedElectronicIndices(N);
    permutedNuclearIndices = allParticlePermutation.indices().head(M);
    permutedElectronicIndices = allParticlePermutation.indices().tail(N).array() - M;

    return {Eigen::PermutationMatrix<Eigen::Dynamic>(permutedNuclearIndices),
            Eigen::PermutationMatrix<Eigen::Dynamic>(permutedElectronicIndices)};
}

namespace YAML {
    Node convert<MolecularGeometry>::encode(const MolecularGeometry &rhs) {
        Node node;
        node["Atoms"] = rhs.atoms();
        node["Electrons"] = rhs.electrons();
        return node;
    }
    bool convert<MolecularGeometry>::decode(const Node &node, MolecularGeometry &rhs) {
        if (!node.IsMap())
            return false;
        rhs = MolecularGeometry({node["Atoms"].as<AtomsVector>(),node["Electrons"].as<ElectronsVector>()});
        return true;
    }

    Emitter &operator<<(Emitter &out, const MolecularGeometry &p) {
        out << BeginMap
            << Key << "Atoms" << Value << p.atoms() << Newline
            << Key << "Electrons" << Value <<  p.electrons()
            << EndMap;
        return out;
    }
}
