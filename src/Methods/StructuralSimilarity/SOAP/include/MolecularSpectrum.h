// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MOLECULARSPECTRUM_H
#define INPSIGHTS_MOLECULARSPECTRUM_H

#include "NeighborhoodExpansion.h"
#include "MolecularGeometry.h"
#include "ParticleKit.h"
#include "NeighborhoodExpander.h"
#include <vector>

namespace SOAP {
    class MolecularSpectrum {
    public:
        MolecularSpectrum() = default;

        explicit MolecularSpectrum(MolecularGeometry molecule);

        void flipSpins();

        MolecularGeometry molecule_;
        MolecularCenters molecularCenters_;
    };
}

#endif //INPSIGHTS_MOLECULARSPECTRUM_H
