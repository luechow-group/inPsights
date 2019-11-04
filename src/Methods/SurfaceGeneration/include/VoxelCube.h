/* Copyright (C) 2019 Michael Heuer.
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

#ifndef INPSIGHTS_VOXELCUBE_H
#define INPSIGHTS_VOXELCUBE_H

#include <vector>
#include <NaturalConstants.h>
#include <Eigen/Core>
#include <DualMC.h>
#include <cstdint>

class VoxelCube {
public:
    using VolumeDataType = uint32_t ;
    using IndexType = dualmc::QuadIndexType;
    using VertexComponentsType = dualmc::VertexComponentsType;

    explicit VoxelCube(
            IndexType dimension = 16,
            VertexComponentsType length = VertexComponentsType(8 * ConversionFactors::angstrom2bohr),
            const Eigen::Matrix<VertexComponentsType,3,1>& origin = {0,0,0},
            bool boxSmoothQ = false);

    std::size_t index(IndexType i, IndexType j, IndexType k) const;

    Eigen::Vector3i getVoxelIndices(const Eigen::Vector3d& pos);

    void add(const Eigen::Vector3d& pos, VolumeDataType weight = 1);

    void add(IndexType i, IndexType j, IndexType k, VolumeDataType weight = 1);

    void shiftDualMCResults(std::vector<dualmc::Vertex>& vertices);

    IndexType getDimension() const;

    VertexComponentsType getLength() const;

    const Eigen::Matrix<VertexComponentsType, 3, 1> &getOrigin() const;

    const std::vector<VolumeDataType> &getData() const;

    VolumeDataType getData(IndexType i, IndexType j, IndexType k) const;

    void smooth(IndexType neighbors);

    VolumeDataType cubeAverage(IndexType i, IndexType j, IndexType k, IndexType neighbors);

    void setData(const std::vector<VolumeDataType> &data);
   
    static constexpr VertexComponentsType offset_ = 0.5;

    bool smoothQ_;
    IndexType dimension_;
    VolumeDataType insideWeight_, totalWeight_;
    VertexComponentsType length_, halfLength_, inverseDimension_;
    Eigen::Matrix<VertexComponentsType,3,1> origin_;
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
