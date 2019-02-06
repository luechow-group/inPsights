//
// Created by Michael Heuer on 2019-02-02.
//

#ifndef INPSIGHTS_MOTIFANALYSIS_H
#define INPSIGHTS_MOTIFANALYSIS_H

#include <vector>
#include "ParticlesVector.h"
#include <GraphAnalysis.h>
#include <list>
#include <Metrics.h>
#include <MolecularGeometry.h>
#include <SpinPairClassification.h>
#include <Varname.h>
#include <spdlog/spdlog.h>

namespace MotifAnalysis{
    enum MotifType {unassigned=0, Core, Valence};

   std::string toString(MotifType type);
   MotifType fromString(const std::string& string);


    class Motif{
    public:
        Motif(const std::list<Eigen::Index>& electronIndices, MotifType type = MotifType::unassigned);

        bool containsQ(Eigen::Index i) const;

        // needed for maps
        bool operator<(const Motif &rhs) const;
        bool operator>(const Motif &rhs) const;
        bool operator<=(const Motif &rhs) const;
        bool operator>=(const Motif &rhs) const;

        MotifType type() const;
        void setType(MotifType type_);

        const std::list<Eigen::Index> &electronIndices() const;
        void setElectronIndices(const std::list<Eigen::Index> &electronIndices_);

    private:
        MotifType type_;
        std::list<Eigen::Index> electronIndices_;
    };

    bool coreElectronQ(const Electron& e, const AtomsVector& atoms, double threshold = 0.01);


    class Motifs{
    public:
        Motifs(const Eigen::MatrixXb &adjacencyMatrix);
        Motifs(const Eigen::MatrixXb &adjacencyMatrix, const MolecularGeometry& molecule);

        Motifs(std::vector<Motif> motifs);


        void classifyMotifs(const MolecularGeometry& molecule);

        void splitCoreMotifs(const MolecularGeometry& molecule);

        void sort();

        std::vector<Motif> motifVector;

    static std::vector<Motif> motifsFromAdjacencyMatrix(const Eigen::MatrixXb &adjacencyMatrix);
    };

    //TODO classification thresholds into settings
}

namespace YAML {
    class Node; class Emitter;
    template<> struct convert<MotifAnalysis::Motif> {
        static Node encode(const MotifAnalysis::Motif &rhs);
        static bool decode(const Node &node, MotifAnalysis::Motif &rhs);
    };
    Emitter &operator<<(Emitter &out, const MotifAnalysis::Motif &rhs) ;
}

#endif //INPSIGHTS_MOTIFANALYSIS_H
