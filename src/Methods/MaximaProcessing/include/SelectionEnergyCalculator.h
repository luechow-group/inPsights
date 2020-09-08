// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
        SelectionInteractionEnergies() = default;

        EnergyResultsBundle<VectorStatistics> intraEnergies;
        EnergyResultsBundle<TriangularMatrixStatistics> interEnergies;
    };

    SelectionEnergyCalculator(
            const std::vector<Sample> &samples,
            const std::vector<DynamicMolecularSelection> &selections);

    void addTopLevel(const Cluster &cluster);

    std::vector<MolecularSelection> getMolecularSelection(const Cluster &cluster) const;

    void add(const Cluster &cluster);

private:
    const std::vector<Sample> &samples_;
    const std::vector<DynamicMolecularSelection>& dynamicSelections_;
public:
    std::vector<MolecularSelection> molecularSelections_;
    SelectionInteractionEnergies selectionInteractions_;
};

namespace YAML {
    class Emitter;
    Emitter& operator<< (Emitter& out, const SelectionEnergyCalculator::SelectionInteractionEnergies& rhs);
}

#endif //INPSIGHTS_SELECTIONENERGYCALCULATOR_H
