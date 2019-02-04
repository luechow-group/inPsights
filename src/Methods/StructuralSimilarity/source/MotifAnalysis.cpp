//
// Created by Michael Heuer on 2019-02-02.
//

#include "MotifAnalysis.h"


namespace MotifAnalysis {

    Motif::Motif(const std::list<Eigen::Index> &electronIndices)
            : electronIndices(electronIndices) {};

    bool Motif::containsQ(Eigen::Index i) const {
        return std::find(electronIndices.begin(), electronIndices.end(), i) != electronIndices.end();
    }

    bool Motif::operator<(const Motif &rhs) const {
        return electronIndices < rhs.electronIndices;
    }

    bool Motif::operator>(const Motif &rhs) const {
        return rhs < *this;
    }

    bool Motif::operator<=(const Motif &rhs) const {
        return !(rhs < *this);
    }

    bool Motif::operator>=(const Motif &rhs) const {
        return !(*this < rhs);
    }


    Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix)
            : motifVector() {
        auto electronIndicesLists = GraphAnalysis::findGraphClusters(adjacencyMatrix);
        for (const auto &list : electronIndicesLists) {
            motifVector.emplace_back(list);
        }
    };

    Motifs::Motifs(std::vector<Motif> motifs)
            : motifVector(motifs) {};
}


namespace YAML {
    Node convert<MotifAnalysis::Motif>::encode(const MotifAnalysis::Motif &rhs) {
        Node node;
        node = rhs.electronIndices;
        return node;
    }

    bool convert<MotifAnalysis::Motif>::decode(const Node &node, MotifAnalysis::Motif &rhs) {
        if (!node.IsSequence())
            return false;

        rhs = node.as<std::list<Eigen::Index>>();
        return true;
    }

    Emitter &operator<<(Emitter &out, const MotifAnalysis::Motif &rhs) {
        out << YAML::Flow << BeginSeq;
        for(auto i : rhs.electronIndices)
            out <<  i;
        out << EndSeq;
        return out;
    };
}