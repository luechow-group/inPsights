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

#ifndef INPSIGHTS_CLUSTERDATA_H
#define INPSIGHTS_CLUSTERDATA_H

#include <Statistics.h>
#include <ParticlesVector.h>
#include <VoxelCube.h>
#include <Motifs.h>
#include <EnergyStatistics.h>

class ClusterData {
public:
    ClusterData() = default;

    ClusterData(unsigned totalNumberOfStructures,
                const std::vector<ElectronsVector> & exemplaricStructures,
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
                const std::vector<VoxelCube>& voxelCubes
                );

    ElectronsVector representativeStructure() const;

    unsigned N_;
    std::vector<ElectronsVector> exemplaricStructures_;
    Motifs motifs_;
    SingleValueStatistics valueStats_, EtotalStats_;
    VectorStatistics EeStats_, intraMotifEnergyStats_; //TODO remove EeStats
    EnergyStatistics::ElectronicEnergy electronicEnergyStats_;
    TriangularMatrixStatistics SeeStats_, interMotifEnergyStats_, ReeStats_;
    MatrixStatistics RenStats_;
    std::vector<VoxelCube> voxelCubes_;
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
