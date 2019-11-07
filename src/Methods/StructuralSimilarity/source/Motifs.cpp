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
#include <spdlog/spdlog.h>
#include <Motifs.h>
#include <ElementInfo.h>


Motifs::Motifs()
: motifVector_({}) {}

Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix)
        : Motifs(motifsFromAdjacencyMatrix(adjacencyMatrix)) {};

Motifs::Motifs(const std::vector<Motif>& motifs)
        : motifVector_(motifs) {};

Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix, const MolecularGeometry & molecule)
        : Motifs(adjacencyMatrix) {
    classifyMotifs(molecule);
};

std::vector<Motif> Motifs::motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix){
    auto lists = GraphAnalysis::findGraphClusters(adjacencyMatrix);

    std::vector<Motif> motifVector;
    for (const auto &list : lists)
        motifVector.emplace_back(list);

    return motifVector;
}

void Motifs::classifyMotifs(const MolecularGeometry& molecule) {

    std::vector<Motif> newMotifVector{};

    for(auto& motif : motifVector_) {
        std::list<Eigen::Index> involvedCores;

        for (const auto & e : motif.electronIndices()) {
            auto[atCoreQ, coreIndex] = molecule.coreElectronQ(e);

            if (atCoreQ)
                involvedCores.push_back(coreIndex);
        }
        involvedCores.sort();
        involvedCores.unique();

        motif.setAtomIndices(involvedCores);
    }

    for(auto& motif : motifVector_) {

        if(motif.atomIndices().empty()) {
            motif.setType(MotifType::Valence);
            newMotifVector.emplace_back(motif);
        } else {
            for(const auto & k : motif.atomIndices()) {
                auto electronsAtCore = molecule.coreElectronsIndices(k);

                //! assert that electronsAtCore is subset of motif.electronIndices()
                if(electronsAtCore.size() == 2) {
                    if (molecule.atoms()[k].type() == Elements::ElementType::H
                    || molecule.atoms()[k].type() == Elements::ElementType::He) {

                        newMotifVector.emplace_back(Motif({}, {k}, MotifType::Core));
                        newMotifVector.emplace_back(Motif(electronsAtCore, {}, MotifType::Valence));
                    } else {
                        newMotifVector.emplace_back(Motif(electronsAtCore, {k}, MotifType::Core));
                    }
                } else if(electronsAtCore.size() == 1) {
                    if (molecule.atoms()[k].type() == Elements::ElementType::H
                    || molecule.atoms()[k].type() == Elements::ElementType::He) {

                        newMotifVector.emplace_back(Motif({}, {k}, MotifType::Core));
                        newMotifVector.emplace_back(Motif(motif.electronIndices(), {}, MotifType::Valence));
                    } else {
                        spdlog::warn("Unexpected motif: only one electron at non-hydrogen core {}",
                                Elements::ElementInfo::symbol(molecule.atoms()[k].type()));
                    }
                }
            }
        }
    }
    motifVector_ = newMotifVector;
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