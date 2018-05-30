//
// Created by Michael Heuer on 25.05.18.
//

#ifndef AMOLQCPP_MOLECULARSPECTRUM_H
#define AMOLQCPP_MOLECULARSPECTRUM_H

#include "NeighborhoodExpansion.h"
#include "MolecularGeometry.h"
#include "ParticleKit.h"
#include "NeighborhoodExpander.h"
#include <vector>

class MolecularSpectrum{
public:
    explicit MolecularSpectrum(MolecularGeometry molecule);

    MolecularGeometry molecule_;
    MolecularCenters molecularCenters_;
};

#endif //AMOLQCPP_MOLECULARSPECTRUM_H
