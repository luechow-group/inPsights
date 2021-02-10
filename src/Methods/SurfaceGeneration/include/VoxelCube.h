// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_VOXELCUBE_H
#define INPSIGHTS_VOXELCUBE_H

#include <vector>
#include <NaturalConstants.h>
#include <Eigen/Core>
#include <DualMC.h>
#include <cstdint>
#include <string>

class VoxelCube {
public:
    using VolumeDataType = uint64_t ;
    using IndexType = dualmc::QuadIndexType;
    using VertexComponentsType = dualmc::VertexComponentsType;

    explicit VoxelCube(
            IndexType dimension = 16,
            VertexComponentsType length = VertexComponentsType(8 * ConversionFactors::angstrom2bohr),
            const Eigen::Matrix<VertexComponentsType,3,1>& origin = {0,0,0},
            bool boxSmoothQ = false);

    VoxelCube(const Eigen::Matrix<IndexType, 3, 1> &dimensions,
              const Eigen::Matrix<VertexComponentsType, 3, 1> &lengths,
              const Eigen::Matrix<VertexComponentsType, 3, 1> &center,
              bool smoothQ);

    std::size_t index(IndexType i, IndexType j, IndexType k) const;

    Eigen::Vector3i getVoxelIndices(const Eigen::Vector3d& pos);

    void add(const Eigen::Vector3d& pos, VolumeDataType weight = 1);

    void add(IndexType i, IndexType j, IndexType k, VolumeDataType weight = 1);

    void shiftDualMCResults(std::vector<dualmc::Vertex>& vertices);

    Eigen::Matrix<IndexType,3,1> getDimensions() const;

    Eigen::Matrix<VertexComponentsType,3,1> getLengths() const;

    const Eigen::Matrix<VertexComponentsType, 3, 1> &getCenter() const;

    const std::vector<VolumeDataType> &getData() const;

    VolumeDataType getData(IndexType i, IndexType j, IndexType k) const;

    void smooth(IndexType neighbors);

    VolumeDataType cubeAverage(IndexType i, IndexType j, IndexType k, IndexType neighbors);

    void setData(const std::vector<VolumeDataType> &data);

    void exportMacmolplt(const std::string& filename, const std::string& comment);

    void setTotalWeight(const double &totalWeight);
   
    static constexpr VertexComponentsType offset_ = 0.5;

    bool smoothQ_;
    Eigen::Matrix<IndexType,3,1> dimensions_;
    VolumeDataType insideWeight_, totalWeight_;
    Eigen::Matrix<VertexComponentsType,3,1> lengths_, inverseDimensions_;
    Eigen::Matrix<VertexComponentsType,3,1> center_;
    std::vector<VolumeDataType> data_;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<VoxelCube> {
        static Node encode(const VoxelCube &rhs);
        static bool decode(const Node &node, VoxelCube &rhs);
    };
    Emitter &operator<<(Emitter &out, const VoxelCube &p) ;
}

#endif //INPSIGHTS_VOXELCUBE_H
