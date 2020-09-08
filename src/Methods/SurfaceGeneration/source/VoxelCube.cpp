// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <VoxelCube.h>
#include <yaml-cpp/yaml.h>
#include <spdlog/spdlog.h>
#include <EigenYamlConversion.h>
#include <Varname.h>

VoxelCube::VoxelCube(
        IndexType dimension,
        VertexComponentsType length,
        const Eigen::Matrix<VertexComponentsType, 3, 1> &origin,
        bool boxSmoothQ)
        :
        smoothQ_(boxSmoothQ),
        dimension_(dimension),
        insideWeight_(0),
        totalWeight_(0),
        length_(length),
        halfLength_(length / 2.0f),
        inverseDimension_(1.0f / float(dimension - 1)),
        origin_(origin),
        data_(static_cast<size_t>(dimension * dimension * dimension)) {
    assert(dimension > 3 && "The cube dimensions must be greater than 3.");
    assert(length > 0 && "Length must be positive.");
}

std::size_t VoxelCube::index(IndexType i, IndexType j, IndexType k) const {
    return i + dimension_ * (j + dimension_ * k);
}

Eigen::Vector3i VoxelCube::getVoxelIndices(const Eigen::Vector3d& pos){
    Eigen::Vector3i indices;
    indices[0] = IndexType((pos[0] + halfLength_ - origin_[0]) / length_ * dimension_);
    indices[1] = IndexType((pos[1] + halfLength_ - origin_[1]) / length_ * dimension_);
    indices[2] = IndexType((pos[2] + halfLength_ - origin_[2]) / length_ * dimension_);
    return indices;
}

void VoxelCube::add(const Eigen::Vector3d &pos, VolumeDataType weight) {
    auto indices = getVoxelIndices(pos);

    IndexType i = indices[0];
    IndexType j = indices[1];
    IndexType k = indices[2];

    add(i,j,k,weight);
};

void VoxelCube::add(IndexType i, IndexType j, IndexType k, VolumeDataType weight) {
    if(totalWeight_ == std::numeric_limits<VolumeDataType>::max())
        spdlog::critical("Overflow in voxel cell. The maxium value for the {} of {} has been exceeded.",
                VARNAME(VolumeDataType),
                std::numeric_limits<VolumeDataType>::max());

    if ((0 <= i && i < dimension_)
        && (0 <= j && j < dimension_)
        && (0 <= k && k < dimension_)) {
        assert(index(i, j, k) < data_.size() && "The index must fit into the vector");
        data_[index(i, j, k)] += static_cast<VolumeDataType>(weight);
        insideWeight_ += weight;
    }
    totalWeight_ += weight;
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

VoxelCube::VolumeDataType VoxelCube::getData(IndexType i, IndexType j, IndexType k) const{
    return data_[index(i,j,k)];
};

void VoxelCube::setData(const std::vector<VoxelCube::VolumeDataType> &data) {
    VoxelCube::data_ = data;
}


VoxelCube::VolumeDataType
VoxelCube::cubeAverage(IndexType i, IndexType j, IndexType k, IndexType neighbors) {

    VolumeDataType average = 0;

    auto numberOfAveragedBoxes = std::pow(neighbors * 2 + 1, 3);

    for (IndexType l = -neighbors; l <= neighbors; ++l) {
        for (IndexType m = -neighbors; m <= neighbors; ++m) {
            for (IndexType n = -neighbors; n <= neighbors; ++n) {

                IndexType lp = l, mp = m, np = n;

                if (i + l < 0 || i + l >= dimension_) lp = 0;

                if (j + m < 0 || j + m >= dimension_) mp = 0;

                if (k + n < 0 || k + n >= dimension_) np = 0;

                average += data_[index(i + lp, j + mp, k + np)];
            }
        }
    }
    // Attention: division by integer
    average /= numberOfAveragedBoxes; 

    return average;
}

void VoxelCube::smooth(VoxelCube::IndexType neighbors) {

    std::vector<VolumeDataType> dataSmoothed(static_cast<size_t>(dimension_ * dimension_ * dimension_));

    for (IndexType i = 0; i < dimension_; ++i)
        for (IndexType j = 0; j < dimension_; ++j)
            for (IndexType k = 0; k < dimension_; ++k)
                dataSmoothed[index(i, j, k)] = cubeAverage(i, j, k, neighbors);

    data_ = dataSmoothed;
}

namespace YAML {
    Node convert<VoxelCube>::encode(const VoxelCube &rhs) {
        Node node;
        node["smoothed"] = rhs.smoothQ_;
        node["dimension"] = rhs.getDimension();
        node["length"] = rhs.getLength();
        node["origin"] = rhs.getOrigin();
        node["data"] = rhs.getData();
        return node;
    }

    bool convert<VoxelCube>::decode(const Node &node, VoxelCube &rhs) {
        if (!node.IsMap())
            return false;

        bool smoothed = false;
        if(node["smoothed"])
            smoothed = node["smoothed"].as<bool>();
        auto dimension = node["dimension"].as<VoxelCube::IndexType>();
        auto length = node["length"].as<VoxelCube::VertexComponentsType>();
        auto origin = node["origin"].as<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>();
        rhs = VoxelCube(dimension, length, origin, smoothed);
        rhs.setData(node["data"].as<std::vector<VoxelCube::VolumeDataType>>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const VoxelCube &rhs) {
        out << BeginMap
            << Key << "smoothed" << Value << rhs.smoothQ_
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
