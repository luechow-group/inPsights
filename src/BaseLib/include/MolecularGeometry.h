//
// Created by Michael Heuer on 08.05.18.
//

#ifndef AMOLQCPP_MOLECULARGEOMETRY_H
#define AMOLQCPP_MOLECULARGEOMETRY_H

#include "ParticlesVector.h"

class MolecularGeometry{
public:
    MolecularGeometry();
    MolecularGeometry(AtomsVector atoms, ElectronsVector electrons);

    const AtomsVector& atoms() const;
    AtomsVector& atoms();

    const ElectronsVector& electrons() const;
    ElectronsVector & electrons();

    Particle<int> operator[](long i) const;

    std::pair<bool,long> findIndexByNumberedType(const NumberedType<int> &numberedType) const;

    NumberedType<int> findNumberedTypeByIndex(unsigned idx) const;

    long numberOfEntities() const;

private:
    AtomsVector atoms_;
    ElectronsVector electrons_;
};

#endif //AMOLQCPP_MOLECULARGEOMETRY_H
