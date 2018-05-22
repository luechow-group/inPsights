//
// Created by Michael Heuer on 08.05.18.
//

#ifndef AMOLQCPP_MOLECULARGEOMETRY_H
#define AMOLQCPP_MOLECULARGEOMETRY_H

#include <utility>
#include "ParticlesVector.h"

class MolecularGeometry{
public:
    MolecularGeometry();
    MolecularGeometry(AtomsVector atoms, ElectronsVector electrons);

    const AtomsVector& atoms() const;
    AtomsVector& atoms();

    const ElectronsVector& electrons() const;
    ElectronsVector & electrons();

    long numberOfEntities();

private:
    AtomsVector atoms_;
    ElectronsVector electrons_;
};

#endif //AMOLQCPP_MOLECULARGEOMETRY_H
