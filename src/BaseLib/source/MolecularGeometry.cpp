//
// Created by Michael Heuer on 08.05.18.
//

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


NumberedType<int> MolecularGeometry::findNumberedTypeByIndex(unsigned idx) const {
    assert(idx < numberOfEntities() && "The index cannot be greater than the number of particles - 1");

    auto M = atoms().numberOfEntities();
    if(idx < M) {
        return atoms().typesVector().getNumberedTypeByIndex(idx).toIntType();
    } else {
        return electrons().typesVector().getNumberedTypeByIndex(idx-M).toIntType();
    }
}

std::pair<bool,long> MolecularGeometry::findIndexByNumberedType(const NumberedType<int> &numberedType) const {
    if(numberedType.type_ >= int(Spins::first()) && numberedType.type_ <= int(Spins::last())) {
        auto boolIdx = electrons().typesVector().findIndexOfNumberedType(
                NumberedSpin(Spins::spinFromInt(numberedType.type_), numberedType.number_));

        boolIdx.second += atoms().numberOfEntities(); // TODO is this the way it should be?
        return boolIdx;
    } else if(numberedType.type_ >= int(Elements::first()) && numberedType.type_ <= int(Elements::last())) {
        return atoms().typesVector().findIndexOfNumberedType(
                NumberedElement(Elements::elementFromInt(numberedType.type_), numberedType.number_));
    } else {
        return {false,0};
    }
}

bool MolecularGeometry::coreElectronQ(long i, double threshold) const {
    for (Eigen::Index k = 0; k < atoms_.numberOfEntities(); ++k)
        if(Metrics::distance(electrons_[i].position(), atoms_[k].position()) <= threshold)
            return true;

    return false;
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

std::list<long> MolecularGeometry::valenceElectronsIndices(double threshold) const {
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
