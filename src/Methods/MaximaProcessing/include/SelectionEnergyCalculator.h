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

#ifndef INPSIGHTS_SELECTIONENERGYCALCULATOR_H
#define INPSIGHTS_SELECTIONENERGYCALCULATOR_H

#include <Statistics.h>
#include <Sample.h>
#include <Cluster.h>
#include <MolecularSelection.h>
#include <CoulombPotential.h>
#include <EnergyResultsBundle.h>
#include <EnergyPartitioning.h>

class SelectionEnergyCalculator{
public:
    struct SelectionInteractionEnergies {
        EnergyResultsBundle<VectorStatistics> intraEnergies;
        EnergyResultsBundle<TriangularMatrixStatistics> interEnergies;
    };

    SelectionEnergyCalculator(
            const std::vector<Sample> &samples,
            const std::vector<DynamicMolecularSelection> &selections);

    void add(const Cluster &cluster);

private:
    const std::vector<Sample> &samples_;
    const std::vector<DynamicMolecularSelection>& dynamicSelections_;
public:
    SelectionInteractionEnergies selectionInteractions_;
};

namespace YAML {
    class Emitter;

    Emitter& operator<< (Emitter& out, const SelectionEnergyCalculator::SelectionInteractionEnergies& rhs);
}

#endif //INPSIGHTS_SELECTIONENERGYCALCULATOR_H
