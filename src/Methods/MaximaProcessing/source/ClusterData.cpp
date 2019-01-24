//
// Created by Michael Heuer on 22.11.18.
//

#include <ClusterData.h>

ClusterData::ClusterData(unsigned totalNumberOfStructures,
            const std::vector<ElectronsVector>& exemplaricStructures,
            const SingleParticlesStatistics & valueStats,
            const SingleParticlesStatistics & TeStats,
            const SingleParticlesStatistics & EeStats,
            const IntraParticlesStatistics & SeeStats,
            const IntraParticlesStatistics & VeeStats,
            const InterParticlesStatistics & VenStats,
            const std::vector<VoxelCube> &voxelCubes
            )
        :
        N_(totalNumberOfStructures),
        exemplaricStructures_(exemplaricStructures),
        valueStats_(valueStats),
        TeStats_(TeStats),
        EeStats_(EeStats),
        SeeStats_(SeeStats),
        VeeStats_(VeeStats),
        VenStats_(VenStats),
        voxelCubes_(voxelCubes)
        {};

ElectronsVector ClusterData::representativeStructure() const {
    return exemplaricStructures_[0];
}

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
        node["VoxelCubes"] = rhs.voxelCubes_;

        return node;
    }

    bool convert<ClusterData>::decode(const Node &node, ClusterData &rhs) {


        // TODO temporary - remove
        std::vector<VoxelCube> cubes = {};
        if(node["VoxelCubes"])
            if (node["VoxelCubes"].IsSequence() && node["VoxelCubes"].size() > 0)
                cubes = node["VoxelCubes"].as<std::vector<VoxelCube>>();

        rhs = ClusterData(
                node["N"].as<unsigned>(),
                node["Structures"].as<std::vector<ElectronsVector>>(),
                node["ValueRange"].as<SingleParticlesStatistics>(),
                node["Te"].as<SingleParticlesStatistics>(),
                node["Ee"].as<SingleParticlesStatistics>(),
                node["SpinCorrelations"].as<IntraParticlesStatistics>(),
                node["Vee"].as<IntraParticlesStatistics>(),
                node["Ven"].as<InterParticlesStatistics>(),
                        cubes
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
            << Key << "Ee" << Comment("[Eh]") << Value << rhs.EeStats_
            << Key << "Vee" << Comment("[Eh]") << Value << rhs.VeeStats_
            << Key << "Ven" << Comment("[Eh]") << Value << rhs.VenStats_
            << Key << "VoxelCubes" << Value << rhs.voxelCubes_
            << EndMap;

        return out;
    }
}
