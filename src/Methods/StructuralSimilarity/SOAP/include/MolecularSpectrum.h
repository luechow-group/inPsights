//
// Created by Michael Heuer on 25.05.18.
//

#ifndef INPSIGHTS_MOLECULARSPECTRUM_H
#define INPSIGHTS_MOLECULARSPECTRUM_H

#include "NeighborhoodExpansion.h"
#include "MolecularGeometry.h"
#include "ParticleKit.h"
#include "NeighborhoodExpander.h"
#include <vector>

class MolecularSpectrum{
public:
    MolecularSpectrum() = default;
    explicit MolecularSpectrum(MolecularGeometry molecule);

    MolecularGeometry molecule_;
    MolecularCenters molecularCenters_;
};

#endif //INPSIGHTS_MOLECULARSPECTRUM_H
