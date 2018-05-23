//
// Created by Michael Heuer on 08.05.18.
//

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
