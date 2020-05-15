/* Copyright (C) 2019-2020 Michael Heuer.
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

#ifndef INPSIGHTS_MOTIFS_H
#define INPSIGHTS_MOTIFS_H

#include <GraphAnalysis.h>
#include <Motif.h>
#include <MolecularGeometry.h>

class Motifs{
public:

    Motifs();
    Motifs(const Eigen::MatrixXb &adjacencyMatrix);
    Motifs(const Eigen::MatrixXb &adjacencyMatrix,
            const MolecularGeometry& molecule);

    Motifs(std::vector<Motif>  motifs);

    void classifyMotifs(const MolecularGeometry& molecule);

    void mergeMotifs(const std::set<size_t>& indices);

    void sort();

    std::vector<Motif> motifs_;

    static std::vector<Motif> motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix);
};

#endif //INPSIGHTS_MOTIFS_H
