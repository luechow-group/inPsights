// Copyright (C) 2019-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MOTIFS_H
#define INPSIGHTS_MOTIFS_H

#include "GraphAnalysis.h"
#include <Motif.h>
#include <MolecularGeometry.h>

class Motifs {
public:

    Motifs();
    Motifs(const Eigen::MatrixXb &adjacencyMatrix);
    Motifs(const Eigen::MatrixXb &adjacencyMatrix,
            const MolecularGeometry& molecule);

    Motifs(std::vector<Motif>  motifs);

    void classifyMotifs(const MolecularGeometry& molecule);

    void mergeMotifs(const std::set<size_t>& indices);

    std::set<size_t> findMotifMergeIndices(const MolecularGeometry& molecule,
                                           const std::vector<std::vector<size_t>> & nucleiMergeList);

    void sort();

    std::vector<Motif> motifs_;

    static std::vector<Motif> motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix);
};

#endif //INPSIGHTS_MOTIFS_H
