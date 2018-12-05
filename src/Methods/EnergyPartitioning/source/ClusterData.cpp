//
// Created by Michael Heuer on 22.11.18.
//

#include <ClusterData.h>

namespace YAML {
    Node convert<ClusterData>::encode(const ClusterData &rhs) {
        Node node;
        node["N"] = rhs.N_;
        node["ValueRange"] = rhs.valueStats_;
        node["Structures"] = rhs.exemplaricStructures_;
        node["SpinCorrelations"] = rhs.SeeStats_;
        node["Te"] = rhs.TeStats_;
        node["Vee"] = rhs.VeeStats_;
        node["Ven"] = rhs.VenStats_;

        return node;
    }

    bool convert<ClusterData>::decode(const Node &node, ClusterData &rhs) {

        rhs = ClusterData(
                node["N"].as<unsigned>(),
                node["Structures"].as<std::vector<ElectronsVector>>(),
                node["ValueRange"].as<Statistics::RunningStatistics<Eigen::VectorXd,unsigned>>(),
                node["Te"].as<Statistics::RunningStatistics<Eigen::VectorXd,unsigned>>(),
                node["SpinCorrelations"].as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true>>(),
                node["Vee"].as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true>>(),
                node["Ven"].as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned>>()
                );
        
        return true;
    }

    Emitter &operator<<(Emitter &out, const ClusterData &rhs) {
        out << BeginMap
            << Key << "N" << Value << rhs.N_
            << Key << "ValueRange" << Value << Comment("[]") << rhs.valueStats_
            << Key << "Structures" << Comment("[a0]") << Value << rhs.exemplaricStructures_ << Newline
            << Key << "SpinCorrelations" << Comment("[]") << Value << rhs.SeeStats_
            << Key << "Te" << Comment("[Eh]") << Value << rhs.TeStats_
            << Key << "Vee" << Comment("[Eh]") << Value << rhs.VeeStats_
            << Key << "Ven" << Comment("[Eh]") << Value << rhs.VenStats_
            << EndMap;

        return out;
    }
}
