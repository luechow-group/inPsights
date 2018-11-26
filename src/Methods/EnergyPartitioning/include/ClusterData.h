//
// Created by Michael Heuer on 22.11.18.
//

#ifndef INPSIGHTS_CLUSTERDATA_H
#define INPSIGHTS_CLUSTERDATA_H

#include <Statistics.h>
#include <ParticlesVector.h>

class AtomsData{
public:
    AtomsData();

private:
    AtomsVector atomsVector_;
    Eigen::MatrixXd Vnn_;
};

class ClusterData {
public:
    ClusterData(unsigned totalNumberOfStructures,
                std::vector<ElectronsVector> exemplaricStructures,
                Statistics::RunningStatistics<Eigen::VectorXd,unsigned> valueStats,
                Statistics::RunningStatistics<Eigen::VectorXd,unsigned> TeStats,
                Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true> SeeStats,
                Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true> VeeStats,
                Statistics::RunningStatistics<Eigen::MatrixXd,unsigned> VenStats)
    :
    N_(totalNumberOfStructures),
    exemplaricStructures_(std::move(exemplaricStructures)),
    valueStats_(std::move(valueStats)),
    TeStats_(std::move(TeStats)),
    SeeStats_(std::move(SeeStats)),
    VeeStats_(std::move(VeeStats)),
    VenStats_(std::move(VenStats))
    {
    };

    ElectronsVector representativeStructure() const {
        return exemplaricStructures_[0];
    }


    unsigned N_;
    std::vector<ElectronsVector> exemplaricStructures_;
    Statistics::RunningStatistics<Eigen::VectorXd,unsigned> valueStats_,TeStats_;
    Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true> SeeStats_, VeeStats_, VnnStats_;
    Statistics::RunningStatistics<Eigen::MatrixXd,unsigned> VenStats_;
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
