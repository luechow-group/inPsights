//
// Created by Michael Heuer on 2019-02-02.
//

#ifndef INPSIGHTS_MOTIFANALYSIS_H
#define INPSIGHTS_MOTIFANALYSIS_H

#include <vector>
#include "ParticlesVector.h"
#include <GraphAnalysis.h>
#include <list>

namespace MotifAnalysis{


    class Motif{
    public:
        Motif(const std::list<Eigen::Index>& electronIndices);

        bool containsQ(Eigen::Index i) const;

        // needed for maps
        bool operator<(const Motif &rhs) const;
        bool operator>(const Motif &rhs) const;
        bool operator<=(const Motif &rhs) const;
        bool operator>=(const Motif &rhs) const;

        std::list<Eigen::Index> electronIndices;
    };

    class Motifs{
    public:
        Motifs(const Eigen::MatrixXb &adjacencyMatrix);

        Motifs(std::vector<Motif> motifs);

        std::vector<Motif> motifVector;
    };
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
