//
// Created by Michael Heuer on 2019-02-06.
//

#include "Motifs.h"
#include <spdlog/spdlog.h>

Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix)
        : Motifs(motifsFromAdjacencyMatrix(adjacencyMatrix)) {};

Motifs::Motifs(const std::vector<Motif>& motifs)
        : motifVector_(motifs) {};

Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix, const MolecularGeometry & molecule)
        : Motifs(adjacencyMatrix) {
    classifyMotifs(molecule);
    splitCoreMotifs(molecule);
};

std::vector<Motif> Motifs::motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix){
    auto lists = GraphAnalysis::findGraphClusters(adjacencyMatrix);

    std::vector<Motif> motifVector;
    for (const auto &list : lists)
        motifVector.emplace_back(list);

    return motifVector;
}

void Motifs::classifyMotifs(const MolecularGeometry& molecule) {
    for(auto& motif : motifVector_) {
        std::vector<bool> atCore;

        for(auto i : motif.electronIndices())
            atCore.emplace_back(molecule.coreElectronQ(i));

        if(std::all_of(atCore.begin(), atCore.end(), [](bool b){ return b; }))
            motif.setType(MotifType::Core);
        else if(std::none_of(atCore.begin(), atCore.end(), [](bool b){ return b; }))
            motif.setType(MotifType::Valence);
        else {
            YAML::Emitter out; out << molecule.electrons() <<  motif;
            spdlog::warn("Found a motif that is not clearly separable into Valence and Core. "
                         "ElectronsVector:\n{0} Motif:\n{1}", out.c_str());
        }
    }
}

void Motifs::splitCoreMotifs(const MolecularGeometry& molecule) {
    std::vector<Motif> newMotifVector{};

    for(const auto& m : motifVector_) {
        if(m.type() == MotifType::Core){
            for (Eigen::Index k = 0; k < molecule.atoms().numberOfEntities(); ++k) {
                auto atCoreK = molecule.coreElectronsIndices(k);

                std::list<Eigen::Index> intersection;
                std::set_intersection(
                        m.electronIndices().begin(), m.electronIndices().end(),
                        atCoreK.begin(), atCoreK.end(), std::back_inserter(intersection));
                assert(intersection.size() <= 2
                       && "The number of electrons at a nuclear cusp must be less than or equal to 2 due to the antisymmetry principle.");

                newMotifVector.emplace_back(Motif(intersection, {k}, MotifType::Core));
            }
        } else
            newMotifVector.emplace_back(m);
    }
    motifVector_ = newMotifVector;
}

void Motifs::sort(){
    std::sort(std::begin(motifVector_), std::end(motifVector_),
              [] (const auto& lhs, const auto& rhs) {
                  return lhs.type() < rhs.type();
              });
}
