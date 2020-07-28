/* Copyright (C) 2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
