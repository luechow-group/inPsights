//
// Created by Michael Heuer on 08.05.18.
//

#include <utility>
#include "MolecularGeometry.h"

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
