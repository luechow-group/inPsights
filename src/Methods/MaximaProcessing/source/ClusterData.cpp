//
// Created by Michael Heuer on 22.11.18.
//

#include <ClusterData.h>

ClusterData::ClusterData(unsigned totalNumberOfStructures,
            const std::vector<ElectronsVector>& exemplaricStructures,
            const SingleValueStatistics & valueStats,
            const VectorStatistics & TeStats,
            const VectorStatistics & EeStats,
            const TriangularMatrixStatistics & SeeStats,
            const TriangularMatrixStatistics & VeeStats,
            const MatrixStatistics & VenStats,
            const Motifs& motifs,
            const SingleValueStatistics & EtotalStats,
            const VectorStatistics & intraMotifEnergyStats,
            const TriangularMatrixStatistics & interMotifEnergyStats,
            const TriangularMatrixStatistics & ReeStats,
            const MatrixStatistics RenStats,
            const std::vector<VoxelCube> &voxelCubes
            )
        :
        N_(totalNumberOfStructures),
        exemplaricStructures_(exemplaricStructures),
        motifs_(motifs),
        valueStats_(valueStats),
        EtotalStats_(EtotalStats),
        EeStats_(EeStats),
        intraMotifEnergyStats_(intraMotifEnergyStats),
        electronicEnergyStats_(TeStats, VeeStats, VenStats),
        SeeStats_(SeeStats),
        interMotifEnergyStats_(interMotifEnergyStats),
        ReeStats_(ReeStats),
        RenStats_(RenStats),
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
        node["Motifs"] = rhs.motifs_.motifVector_;
        node["Etotal"] = rhs.EtotalStats_;
        node["IntraMotifEnergies"] = rhs.intraMotifEnergyStats_;
        node["InterMotifEnergies"] = rhs.interMotifEnergyStats_;
        node["Structures"] = rhs.exemplaricStructures_;
        node["SpinCorrelations"] = rhs.SeeStats_;
        node["Ree"] = rhs.ReeStats_;
        node["Ren"] = rhs.RenStats_;
        node["Te"] = rhs.electronicEnergyStats_.Te();
        node["Vee"] = rhs.electronicEnergyStats_.Vee();
        node["Ven"] = rhs.electronicEnergyStats_.Ven();
        node["VoxelCubes"] = rhs.voxelCubes_;

        return node;
    }

    bool convert<ClusterData>::decode(const Node &node, ClusterData &rhs) {

        // TODO temporary - remove
        std::vector<VoxelCube> cubes = {};
        if(node["VoxelCubes"])
            if (node["VoxelCubes"].IsSequence() && node["VoxelCubes"].size() > 0)
                cubes = node["VoxelCubes"].as<std::vector<VoxelCube>>();

        auto motifVector = node["Motifs"].as<std::vector<Motif>>();

        rhs = ClusterData(
                node["N"].as<unsigned>(),
                node["Structures"].as<std::vector<ElectronsVector>>(),
                node["ValueRange"].as<SingleValueStatistics>(),
                node["Te"].as<VectorStatistics>(),
                node["Ee"].as<VectorStatistics>(),
                node["SpinCorrelations"].as<TriangularMatrixStatistics>(),
                node["Vee"].as<TriangularMatrixStatistics>(),
                node["Ven"].as<MatrixStatistics>(),
                Motifs(motifVector),
                node["Etotal"].as<SingleValueStatistics>(),
                node["IntraMotifEnergies"].as<VectorStatistics>(),
                node["InterMotifEnergies"].as<TriangularMatrixStatistics>(),
                node["Ree"].as<TriangularMatrixStatistics>(),
                node["Ren"].as<MatrixStatistics>(),
                cubes
                );
        
        return true;
    }

    Emitter &operator<<(Emitter &out, const ClusterData &rhs) {
        out << BeginMap
            << Key << "N" << Value << rhs.N_
            << Key << "ValueRange" << Value << Comment("[]") << rhs.valueStats_
            << Key << "Motifs" << Value << rhs.motifs_.motifVector_
            << Key << "Etotal" << Comment("[Eh]") << Value << rhs.EtotalStats_
            << Key << "IntraMotifEnergies" << Comment("[Eh]") << Value << rhs.intraMotifEnergyStats_
            << Key << "InterMotifEnergies" << Comment("[Eh]") << Value << rhs.interMotifEnergyStats_
            << Key << "Structures" << Comment("[a0]") << Value << rhs.exemplaricStructures_ << Newline
            << Key << "SpinCorrelations" << Comment("[]") << Value << rhs.SeeStats_
            << Key << "Ree" << Comment("[a0]") << Value << rhs.ReeStats_
            << Key << "Ren" << Comment("[a0]") << Value << rhs.RenStats_
            << Key << "Te" << Comment("[Eh]") << Value << rhs.electronicEnergyStats_.Te()
            << Key << "Ee" << Comment("[Eh]") << Value << rhs.EeStats_
            << Key << "Vee" << Comment("[Eh]") << Value << rhs.electronicEnergyStats_.Vee()
            << Key << "Ven" << Comment("[Eh]") << Value << rhs.electronicEnergyStats_.Ven()
            << Key << "VoxelCubes" << Value << rhs.voxelCubes_
            << EndMap;

        return out;
    }
}
