// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CLUSTERDATA_H
#define INPSIGHTS_CLUSTERDATA_H

#include <Statistics.h>
#include <ParticlesVector.h>
#include <VoxelCube.h>
#include <Motifs.h>
#include <EnergyStatistics.h>
#include <SelectionEnergyCalculator.h>

class ClusterData {
public:
    ClusterData() = default;

    ClusterData(unsigned totalNumberOfStructures,
                const std::vector<ElectronsVector> & exemplaricStructures,
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
                const std::vector<VoxelCube>& seds,
                const Eigen::MatrixXd& sedOverlaps
    );

    // TODO Refactor (hacky solution for local clustering)
    ClusterData(unsigned totalNumberOfStructures,
                const std::vector<ElectronsVector> & exemplaricStructures,
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
                const std::vector<VoxelCube>& seds,
                const Eigen::MatrixXd& sedOverlaps,
                const SelectionEnergyCalculator::SelectionInteractionEnergies & selectionInteractionEnergies,
                const std::vector<MolecularSelection>& selections
                );

    ElectronsVector representativeStructure() const;

    unsigned N_;
    std::vector<ElectronsVector> exemplaricStructures_;
    ElectronsVector sampleAverage_;
    Motifs motifs_;
    SingleValueStatistics valueStats_, EtotalStats_;
    VectorStatistics EeStats_, intraMotifEnergyStats_; //TODO remove EeStats
    EnergyStatistics::ElectronicEnergy electronicEnergyStats_;
    TriangularMatrixStatistics SeeStats_, interMotifEnergyStats_, ReeStats_;
    MatrixStatistics RenStats_;
    std::vector<VoxelCube> voxelCubes_;
    Eigen::MatrixXd overlaps_;
    SelectionEnergyCalculator::SelectionInteractionEnergies  selectionInteractionEnergies_;
    std::vector<MolecularSelection> selections_;
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
