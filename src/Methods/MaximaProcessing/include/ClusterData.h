//
// Created by Michael Heuer on 22.11.18.
//

#ifndef INPSIGHTS_CLUSTERDATA_H
#define INPSIGHTS_CLUSTERDATA_H

#include <Statistics.h>
#include <ParticlesVector.h>

class ClusterData {
public:
    ClusterData() = default;

    ClusterData(unsigned totalNumberOfStructures,
                std::vector<ElectronsVector> exemplaricStructures,
                SingleParticlesStatistics valueStats,
                SingleParticlesStatistics TeStats,
                SingleParticlesStatistics EeStats,
                IntraParticlesStatistics SeeStats,
                IntraParticlesStatistics VeeStats,
                InterParticlesStatistics VenStats);

    ElectronsVector representativeStructure() const;

    unsigned N_;
    std::vector<ElectronsVector> exemplaricStructures_;
    SingleParticlesStatistics valueStats_,TeStats_, EeStats_;
    IntraParticlesStatistics SeeStats_, VeeStats_;
    InterParticlesStatistics VenStats_;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<ClusterData> {
        static Node encode(const ClusterData &rhs);
        static bool decode(const Node &node, ClusterData &rhs);
    };
    Emitter &operator<<(Emitter &out, const ClusterData &p) ;
}

#endif //INPSIGHTS_CLUSTERDATA_H
