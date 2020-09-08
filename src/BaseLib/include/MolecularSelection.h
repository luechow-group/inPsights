// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MOLECULARSELECTION_H
#define INPSIGHTS_MOLECULARSELECTION_H

#include "MolecularSelection.h"
#include <ParticleIndices.h>
#include <ParticleSelection.h>

struct MolecularSelection {
    MolecularSelection() = default;
    MolecularSelection(ParticleIndices electronIndices);
    MolecularSelection(ParticleIndices electronIndices, ParticleIndices nucleiIndices);


    bool operator==(const MolecularSelection& rhs);

    ParticleIndices electrons_;
    ParticleIndices nuclei_;
};

namespace YAML {
    class Emitter;

    Emitter& operator<< (Emitter& out, const MolecularSelection& rhs);
}

struct DynamicMolecularSelection {
    DynamicMolecularSelection(Settings::ParticleSelection particleSelectionSettings, ParticleIndices nucleiIndices);

    ParticleIndices selectNearestElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei) const;

    MolecularSelection toMolecularSelection(const ElectronsVector &electrons, const AtomsVector &nuclei) const;

    Settings::ParticleSelection particleSelectionSettings_;
    ParticleIndices nuclei_;
};

#endif //INPSIGHTS_MOLECULARSELECTION_H
