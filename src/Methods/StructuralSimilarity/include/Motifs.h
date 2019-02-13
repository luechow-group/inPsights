//
// Created by Michael Heuer on 2019-02-06.
//

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

    Motifs(const std::vector<Motif>& motifs);


    void classifyMotifs(const MolecularGeometry& molecule);

    void splitCoreMotifs(const MolecularGeometry& molecule);

    void sort();

    std::vector<Motif> motifVector_;

    static std::vector<Motif> motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix);
};

#endif //INPSIGHTS_MOTIFS_H
