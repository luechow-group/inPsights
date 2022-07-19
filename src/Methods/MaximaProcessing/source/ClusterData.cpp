// Copyright (C) 2018-2019 Michael Heuer.
// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

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
                         const Eigen::MatrixXd& sedOverlaps,
                         const std::vector<unsigned>& subCounts,
                         const std::vector<SingleValueStatistics>& subValueStats,
                         const std::vector<Eigen::VectorXd> &eigenvalues,
                         const std::vector<std::vector<Eigen::Vector3d>> &eigenvectors
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
        overlaps_(sedOverlaps),
        subN_(subCounts),
        subValueStats_(subValueStats),
        eigenvalues_(eigenvalues),
        eigenvectors_(eigenvectors)
{};

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
            const Eigen::MatrixXd& sedOverlaps,
            const SelectionEnergyCalculator::SelectionInteractionEnergies & selectionInteractionEnergies,
            const std::vector<MolecularSelection>& selections,
            const std::vector<unsigned>& subCounts,
            const std::vector<SingleValueStatistics>& subValueStats
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
        overlaps_(sedOverlaps),
        selectionInteractionEnergies_(selectionInteractionEnergies),
        selections_(selections),
        subN_(subCounts),
        subValueStats_(subValueStats),
        eigenvalues_(),
        eigenvectors_()
        {};

ElectronsVector ClusterData::representativeStructure() const {
    return exemplaricStructures_[0];
}

namespace YAML {
    Node convert<ClusterData>::encode(const ClusterData &rhs) {
        Node node;
        node["N"] = rhs.N_;
        node["ValueRange"] = rhs.valueStats_;
        node["Motifs"] = rhs.motifs_.motifs_;
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
        node["SubStructureN"] = rhs.subN_;
        node["SubStructureValueRange"] = rhs.subValueStats_;
        return node;
    }

    bool convert<ClusterData>::decode(const Node &node, ClusterData &rhs) {

        std::vector<VoxelCube> cubes = {};
        if(node["VoxelCubes"])
            if (node["VoxelCubes"].IsSequence() && node["VoxelCubes"].size() > 0)
                cubes = node["VoxelCubes"].as<std::vector<VoxelCube>>();

        Eigen::MatrixXd sedOverlaps;
        if(node["SedOverlaps"])
            if (node["SedOverlaps"].IsMap()  && node["SedOverlaps"][0])
                sedOverlaps= node["SedOverlaps"].as<Eigen::MatrixXd>();

        std::vector<Motif> motifVector;
        if (node["Motifs"])
            motifVector = node["Motifs"].as<std::vector<Motif>>();

        auto structures = node["Structures"].as<std::vector<ElectronsVector>>();

        ElectronsVector sampleAverage;
        if(node["SampleAverage"])
            sampleAverage = node["SampleAverage"].as<ElectronsVector>();
        else
            sampleAverage = {};

        std::vector<unsigned> subCounts;
        if (node["SubStructureN"]) {
            subCounts = node["SubStructureN"].as<std::vector<unsigned>>();
        }
        else{
            for (int i=0; i<structures.size();++i){
                subCounts.emplace_back(NAN);
            }
        }
        std::vector<SingleValueStatistics> subValueStats;
        if (node["SubStructureValueRange"]) {
            subValueStats = node["SubStructureValueRange"].as<std::vector<SingleValueStatistics>>();
        }

        std::vector<Eigen::VectorXd> eigenvalues;
        std::vector<std::vector<Eigen::Vector3d>> eigenvectors;
        for (auto subNode : node["Structures"]){
            if (subNode["Eigenvalues"]){
                eigenvalues.emplace_back(subNode["Eigenvalues"].as<Eigen::VectorXd>());
                std::vector<Eigen::Vector3d> vecs3d;
                for (auto eigenVec : subNode["Eigenvectors"]){
                    vecs3d.emplace_back(eigenVec.as<Eigen::Vector3d>());
                }
                eigenvectors.emplace_back(vecs3d);
            }
        }

        VectorStatistics TeStats;
        if (node["Te"])
            TeStats = node["Te"].as<VectorStatistics>();
        VectorStatistics EeStats;
        if (node["Ee"])
            EeStats = node["Te"].as<VectorStatistics>();
        TriangularMatrixStatistics SpinCorrelationsStats;
        if (node["SpinCorrelations"])
            SpinCorrelationsStats = node["SpinCorrelations"].as<TriangularMatrixStatistics>();
        TriangularMatrixStatistics VeeStats;
        if (node["Vee"])
            VeeStats = node["Vee"].as<TriangularMatrixStatistics>();
        MatrixStatistics VenStats;
        if (node["Ven"])
            VenStats = node["Ven"].as<MatrixStatistics>();
        SingleValueStatistics EtotalStats;
        if (node["Etotal"])
            EtotalStats = node["Ven"].as<SingleValueStatistics>();
        VectorStatistics IntraMotifEnergiesStats;
        if (node["IntraMotifEnergies"])
            IntraMotifEnergiesStats = node["IntraMotifEnergies"].as<VectorStatistics>();
        TriangularMatrixStatistics InterMotifEnergiesStats;
        if (node["InterMotifEnergies"])
            InterMotifEnergiesStats = node["InterMotifEnergies"].as<TriangularMatrixStatistics>();
        TriangularMatrixStatistics ReeStats;
        if (node["Ree"])
            ReeStats = node["Ree"].as<TriangularMatrixStatistics>();
        MatrixStatistics RenStats;
        if (node["Ren"])
            RenStats = node["Ren"].as<MatrixStatistics>();

        rhs = ClusterData(
                node["N"].as<unsigned>(),
                structures,
                sampleAverage,
                node["ValueRange"].as<SingleValueStatistics>(),
                TeStats,
                EeStats,
                SpinCorrelationsStats,
                VeeStats,
                VenStats,
                Motifs(motifVector),
                EtotalStats,
                IntraMotifEnergiesStats,
                InterMotifEnergiesStats,
                ReeStats,
                RenStats,
                cubes,
                sedOverlaps,
                subCounts,
                subValueStats,
                eigenvalues,
                eigenvectors
                );
        
        return true;
    }

    Emitter &operator<<(Emitter &out, const ClusterData &rhs) {
        out << BeginMap
            << Key << "N" << Value << rhs.N_
            << Key << "ValueRange" << Value << Comment("[]") << rhs.valueStats_
            << Key << "SampleAverage" << Comment("[a0]") << Value << rhs.sampleAverage_<< Newline;
        if (rhs.EtotalStats_.getTotalWeight() > 0)
            out << Key << "Etotal" << Comment("[Eh]") << Value << rhs.EtotalStats_;
        if (not rhs.motifs_.motifs_.empty())
            out << Key << "Motifs" << Value << rhs.motifs_.motifs_
                << Key << "IntraMotifEnergies" << Comment("[Eh]") << Value << rhs.intraMotifEnergyStats_
                << Key << "InterMotifEnergies" << Comment("[Eh]") << Value << rhs.interMotifEnergyStats_;
        if (not rhs.selections_.empty())
            out << Key << "Selections" << Value << rhs.selections_
                << Key << "SelectionEnergyCalculation" << Value << rhs.selectionInteractionEnergies_;
        out << Key << "Structures" << Comment("[a0]") << Value << rhs.exemplaricStructures_ << Newline;
        if (rhs.SeeStats_.getTotalWeight() > 0)
            out << Key << "SpinCorrelations" << Comment("[]") << Value << rhs.SeeStats_;
        if (rhs.ReeStats_.getTotalWeight() > 0)
            out << Key << "Ree" << Comment("[a0]") << Value << rhs.ReeStats_
                << Key << "Ren" << Comment("[a0]") << Value << rhs.RenStats_
                << Key << "Te" << Comment("[Eh]") << Value << rhs.electronicEnergyStats_.Te()
                << Key << "Ee" << Comment("[Eh]") << Value << rhs.EeStats_
                << Key << "Vee" << Comment("[Eh]") << Value << rhs.electronicEnergyStats_.Vee()
                << Key << "Ven" << Comment("[Eh]") << Value << rhs.electronicEnergyStats_.Ven();
        if (not rhs.voxelCubes_.empty())
            out << Key << "VoxelCubes" << Value << rhs.voxelCubes_;
        if (rhs.overlaps_.size() > 0)
            out << Key << "SedOverlaps" << Value << rhs.overlaps_;

        out << Key << "SubStructureN" << Value << rhs.subN_;

        out << Key << "SubStructureValueRange" << Value << rhs.subValueStats_
            << EndMap;

        return out;
    }
}
