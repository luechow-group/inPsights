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

#include <VoxelCube.h>
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

VoxelCube::VoxelCube(
        IndexType dimension,
        VertexComponentsType length,
        const Eigen::Matrix<VertexComponentsType, 3, 1>& origin)
        :
        dimension_(dimension),
        length_(length),
        halfLength_(length / 2.0f),
        inverseDimension_(1.0f / (dimension - 1)),
        origin_(origin),
        data_(static_cast<unsigned long>(dimension * dimension * dimension)) {
    assert(dimension > 3 && "The cube dimensions must be greater than 3.");
    assert(length > 0 && "Length must be positive.");
}

long VoxelCube::index(IndexType i, IndexType j, IndexType k) {
    return i + dimension_ * (j + dimension_ * k);
}

void VoxelCube::add(const Eigen::Vector3d &pos, IndexType weight) {
    auto i = IndexType((pos[0] + halfLength_ - origin_[0]) / length_ * dimension_);
    auto j = IndexType((pos[1] + halfLength_ - origin_[1]) / length_ * dimension_);
    auto k = IndexType((pos[2] + halfLength_ - origin_[2]) / length_ * dimension_);

    if((0 <= i && i < dimension_)
    && (0 <= j && j < dimension_)
    && (0 <= k && k < dimension_)) {
        assert(index(i, j, k) < data_.size() && "The index must fit into the vector");
        data_[index(i, j, k)] += static_cast<VolumeDataType>(weight);
    }
};

void VoxelCube::shiftDualMCResults(std::vector<dualmc::Vertex> &vertices) {
    for (auto &v : vertices) {
        v.x = (v.x * inverseDimension_ - offset_) * length_ + origin_[0];
        v.y = (v.y * inverseDimension_ - offset_) * length_ + origin_[1];
        v.z = (v.z * inverseDimension_ - offset_) * length_ + origin_[2];
    }
}

VoxelCube::IndexType VoxelCube::getDimension() const {
    return dimension_;
}

VoxelCube::VertexComponentsType VoxelCube::getLength() const {
    return length_;
}

const Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1> &VoxelCube::getOrigin() const {
    return origin_;
}

const std::vector<VoxelCube::VolumeDataType> &VoxelCube::getData() const {
    return data_;
}

void VoxelCube::setData(const std::vector<VoxelCube::VolumeDataType> &data) {
    VoxelCube::data_ = data;
}

namespace YAML {
    Node convert<VoxelCube>::encode(const VoxelCube &rhs) {
        Node node;
        node["dimension"] = rhs.getDimension();
        node["length"] = rhs.getLength();
        node["origin"] = rhs.getOrigin();
        node["data"] = rhs.getData();
        return node;
    }

    bool convert<VoxelCube>::decode(const Node &node, VoxelCube &rhs) {
        if (!node.IsMap())
            return false;

        auto dimension = node["dimension"].as<VoxelCube::IndexType>();
        auto length = node["length"].as<VoxelCube::VertexComponentsType>();
        auto origin = node["origin"].as<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>();
        rhs = VoxelCube(dimension, length, origin);
        rhs.setData(node["data"].as<std::vector<VoxelCube::VolumeDataType>>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const VoxelCube &rhs) {
        out << BeginMap
            << Key << "dimension" << Value << rhs.getDimension()
            << Key << "length" << Value << rhs.getLength()
            << Key << "origin" << Value << rhs.getOrigin()
            << Key << "data" << Value << Flow << BeginSeq;
        for (const auto &i : rhs.getData()) {
            out << int(i);
        }
        out << EndSeq << EndMap;
        return out;
    }
}
