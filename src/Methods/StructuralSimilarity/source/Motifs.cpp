/* Copyright (C) 2019 Michael Heuer.
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

#include "Motifs.h"
#include <Metrics.h>
#include <ElementInfo.h>
#include <MapUtils.h>
#include <spdlog/spdlog.h>
#include <utility>

Motifs::Motifs()
: motifVector_({}) {}

Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix)
        : Motifs(motifsFromAdjacencyMatrix(adjacencyMatrix)) {};

Motifs::Motifs(std::vector<Motif> motifs)
        : motifVector_(std::move(motifs)) {};

Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix, const MolecularGeometry & molecule)
        : Motifs(adjacencyMatrix) {
    classifyMotifs(molecule);
};

std::vector<Motif> Motifs::motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix){
    auto sets = GraphAnalysis::findGraphClusters(adjacencyMatrix);

    std::vector<Motif> motifVector;
    for (const auto &set : sets)
        motifVector.emplace_back(set);

    return motifVector;
}

void Motifs::classifyMotifs(const MolecularGeometry& molecule) {
    std::set<Eigen::Index> remainingNuclei;
    for (long i = 0; i < molecule.atoms().numberOfEntities(); ++i)
        remainingNuclei.emplace(i);

    std::vector<Motif> newMotifsVector;

    // check all motifs
    for(auto& motif : motifVector_) {
        std::set<Eigen::Index> involvedNuclei;

        std::map<long, long> electronToNucleusMap;

        // check for each electron if it is at an nucleus
        for (const auto &electronIndex : motif.electronIndices()) {
            auto[atNucleusQ, nucleusIndex] = molecule.electronAtNucleusQ(electronIndex);

            electronToNucleusMap.emplace(electronIndex, nucleusIndex);

            // electrons at H or He nuclei are still valence electrons so the nucleus is not involved in the motif
            if (atNucleusQ &&
            !(molecule.atoms()[nucleusIndex].type() == Elements::ElementType::H
            || molecule.atoms()[nucleusIndex].type() == Elements::ElementType::He))
                involvedNuclei.emplace(nucleusIndex);
        }

        if(involvedNuclei.empty()) {
            motif.setType(MotifType::Valence);
            newMotifsVector.emplace_back(motif);
        } else {
            // create a core motifs fro every involved nuclei with all electrons in it.
            for(const auto & nucleusIndex : involvedNuclei) {
                auto electronKeys = MapUtils::findByValue(electronToNucleusMap, nucleusIndex);
                assert(!electronKeys.empty() && "Some electrons must be at involved nuclei at this point.");

                newMotifsVector.emplace_back(
                        Motif{std::set<Eigen::Index>(
                                std::begin(electronKeys),
                                std::end(electronKeys)),
                                {nucleusIndex}, MotifType::Core});
            }
        }

        std::set<long> diff{};
        std::set_difference(remainingNuclei.begin(), remainingNuclei.end(), involvedNuclei.begin(),
                            involvedNuclei.end(),
                            std::inserter(diff, diff.begin()));
        remainingNuclei = diff;
    }

    // remaing nuclei (only H, He) being not involved in any motifs are separate core motifs
    for(auto & nucleusIndex : remainingNuclei)
        newMotifsVector.emplace_back(Motif({}, {nucleusIndex}, MotifType::Core));

    motifVector_ = newMotifsVector;

    sort();
}

void Motifs::sort(){
    std::sort(std::begin(motifVector_), std::end(motifVector_),
              [] (const auto& lhs, const auto& rhs) {
        if(lhs.type() < rhs.type())
            return true;
        else if (lhs.type() > rhs.type())
            return false;
        else
            return lhs.electronIndices().size() > rhs.electronIndices().size();
    });
}
