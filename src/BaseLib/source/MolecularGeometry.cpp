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
        return {false,0};
    }
}

std::tuple<bool,Eigen::Index> MolecularGeometry::electronAtNucleusQ(long i, double threshold) const {
    for (Eigen::Index k = 0; k < atoms_.numberOfEntities(); ++k)
        if(Metrics::distance(electrons_[i].position(), atoms_[k].position()) <= threshold)
            return {true, k};

    return {false, 0};
}

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
