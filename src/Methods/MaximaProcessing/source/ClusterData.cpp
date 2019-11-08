/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <ClusterData.h>


ClusterData::ClusterData(unsigned totalNumberOfStructures,
            const std::vector<ElectronsVector>& exemplaricStructures,
            const ElectronsVector& sampleAverage,
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
            const std::vector<VoxelCube> &seds,
            const Eigen::MatrixXd& sedOverlaps
            )
        :
        N_(totalNumberOfStructures),
        exemplaricStructures_(exemplaricStructures),
        sampleAverage_(sampleAverage),
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
        voxelCubes_(seds),
        overlaps_(sedOverlaps)
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
        node["SampleAverage"] = rhs.sampleAverage_;
        node["IntraMotifEnergies"] = rhs.intraMotifEnergyStats_;
        node["InterMotifEnergies"] = rhs.interMotifEnergyStats_;
        node["Structures"] = rhs.exemplaricStructures_;
        node["AverageStructure"] = rhs.exemplaricStructures_;
        node["SpinCorrelations"] = rhs.SeeStats_;
        node["Ree"] = rhs.ReeStats_;
        node["Ren"] = rhs.RenStats_;
        node["Te"] = rhs.electronicEnergyStats_.Te();
        node["Vee"] = rhs.electronicEnergyStats_.Vee();
        node["Ven"] = rhs.electronicEnergyStats_.Ven();
        node["VoxelCubes"] = rhs.voxelCubes_;
        node["SedOverlaps"] = rhs.overlaps_;
        return node;
    }

    bool convert<ClusterData>::decode(const Node &node, ClusterData &rhs) {

        std::vector<VoxelCube> cubes = {};
        if(node["VoxelCubes"])
            if (node["VoxelCubes"].IsSequence() && node["VoxelCubes"].size() > 0)
                cubes = node["VoxelCubes"].as<std::vector<VoxelCube>>();

        Eigen::MatrixXd sedOverlaps;
        if(node["SedOverlaps"])
            if (node["SedOverlaps"].IsMap() && node["SedOverlaps"][0].IsMap())
                sedOverlaps= node["SedOverlaps"].as<Eigen::MatrixXd>();

        auto motifVector = node["Motifs"].as<std::vector<Motif>>();

        auto structures = node["Structures"].as<std::vector<ElectronsVector>>();

        ElectronsVector sampleAverage;
        if(node["SampleAverage"])
            sampleAverage = node["SampleAverage"].as<ElectronsVector>();
        else
            sampleAverage = {};

        rhs = ClusterData(
                node["N"].as<unsigned>(),
                structures,
                sampleAverage,
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
                cubes,
                sedOverlaps
                );
        
        return true;
    }

    Emitter &operator<<(Emitter &out, const ClusterData &rhs) {
        out << BeginMap
            << Key << "N" << Value << rhs.N_
            << Key << "ValueRange" << Value << Comment("[]") << rhs.valueStats_
            << Key << "SampleAverage" << Comment("[a0]") << Value << rhs.sampleAverage_
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
            << Key << "SedOverlaps" << Value << rhs.overlaps_
            << EndMap;

        return out;
    }
}
