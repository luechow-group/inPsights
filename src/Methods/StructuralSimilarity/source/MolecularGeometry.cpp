//
// Created by Michael Heuer on 22.05.18.
//

#include "MolecularGeometry.h"

MolecularGeometry::MolecularGeometry(AtomsVector atoms, ElectronsVector electrons)
        : atoms_(std::move(atoms)),
          electrons_(std::move(electrons))
{};

const AtomsVector& MolecularGeometry::atoms() const { return atoms_; }

AtomsVector& MolecularGeometry::atoms() { return atoms_; }

const ElectronsVector& MolecularGeometry::electrons() const { return electrons_; }

ElectronsVector & MolecularGeometry::electrons() { return electrons_; }

long MolecularGeometry::numberOfEntities(){
    return atoms_.numberOfEntities() + electrons_.numberOfEntities();
}
