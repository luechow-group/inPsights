// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <VoxelCube.h>
#include <yaml-cpp/yaml.h>
#include <spdlog/spdlog.h>
#include <EigenYamlConversion.h>
#include <Varname.h>
#include <Enumerate.h>
#include <ToString.h>
#include <fstream>
#include <cmath>
#include <iostream>

VoxelCube::VoxelCube(const Eigen::Matrix<IndexType, 3, 1> &dimensions,
                     const Eigen::Matrix<VertexComponentsType, 3, 1> &lengths,
                     const Eigen::Matrix<VertexComponentsType, 3, 1> &center,
                     bool smoothQ)
        : smoothQ_(smoothQ), dimensions_(dimensions), insideWeight_(0), totalWeight_(0),
          lengths_(lengths),
          inverseDimensions_({1.0f / float(dimensions[0] - 1),
                              1.0f / float(dimensions[1] - 1),
                              1.0f / float(dimensions[2] - 1)}),
          center_(center),
          data_(static_cast<size_t>(dimensions.prod())) {}

VoxelCube::VoxelCube(
        IndexType dimension,
        VertexComponentsType length,
        const Eigen::Matrix<VertexComponentsType, 3, 1> &center,
        bool boxSmoothQ)
        :
        smoothQ_(boxSmoothQ),
        dimensions_({dimension, dimension, dimension}),
        insideWeight_(0),
        totalWeight_(0),
        lengths_({length, length, length}),
        inverseDimensions_({1.0f / float(dimension - 1), 1.0f / float(dimension - 1), 1.0f / float(dimension - 1)}),
        center_(center),
        data_(static_cast<size_t>(dimension * dimension * dimension)) {
    assert(dimension > 3 && "The cube dimensions must be greater than 3.");
    assert(length > 0 && "Length must be positive.");
}

std::size_t VoxelCube::index(IndexType i, IndexType j, IndexType k) const {
    return i + dimensions_[0] * (j + dimensions_[1] * k);
}

Eigen::Vector3i VoxelCube::getVoxelIndices(const Eigen::Vector3d& pos){
    Eigen::Vector3i indices;
    auto origin = center_ - lengths_/2;
    for (int i=0; i<3; ++i){
        indices[i] = IndexType(lround((pos[i] - origin[i]) / lengths_[i] * dimensions_[i]));
    }
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
        spdlog::critical("Overflow in voxel cell. The maximum value for the {} of {} has been exceeded.",
                VARNAME(VolumeDataType),
                std::numeric_limits<VolumeDataType>::max());

    if ((0 <= i && i < dimensions_[0])
        && (0 <= j && j < dimensions_[1])
        && (0 <= k && k < dimensions_[2])) {
        assert(index(i, j, k) < data_.size() && "The index must fit into the vector");
        data_[index(i, j, k)] += static_cast<VolumeDataType>(weight);
        insideWeight_ += weight;
    }
    totalWeight_ += weight;
};


void VoxelCube::shiftDualMCResults(std::vector<dualmc::Vertex> &vertices) {
    for (auto &v : vertices) {
        v.x = (v.x * inverseDimensions_[0] - offset_) * lengths_[0] + center_[0];
        v.y = (v.y * inverseDimensions_[1] - offset_) * lengths_[1] + center_[1];
        v.z = (v.z * inverseDimensions_[2] - offset_) * lengths_[2] + center_[2];
    }
}

Eigen::Matrix<VoxelCube::IndexType,3,1> VoxelCube::getDimensions() const {
    return dimensions_;
}

Eigen::Matrix<VoxelCube::VertexComponentsType,3,1> VoxelCube::getLengths() const {
    return lengths_;
}

const Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1> &VoxelCube::getCenter() const {
    return center_;
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

                if (i + l < 0 || i + l >= dimensions_[0]) lp = 0;

                if (j + m < 0 || j + m >= dimensions_[1]) mp = 0;

                if (k + n < 0 || k + n >= dimensions_[2]) np = 0;

                average += data_[index(i + lp, j + mp, k + np)];
            }
        }
    }
    // Attention: division by integer
    average /= numberOfAveragedBoxes; 

    return average;
}

void VoxelCube::smooth(VoxelCube::IndexType neighbors) {

    std::vector<VolumeDataType> dataSmoothed(static_cast<size_t>(dimensions_[0] * dimensions_[1] * dimensions_[2]));

    for (IndexType i = 0; i < dimensions_[0]; ++i)
        for (IndexType j = 0; j < dimensions_[1]; ++j)
            for (IndexType k = 0; k < dimensions_[2]; ++k)
                dataSmoothed[index(i, j, k)] = cubeAverage(i, j, k, neighbors);

    data_ = dataSmoothed;
}

void VoxelCube::exportMacmolplt(const std::string& filename, const std::string& comment) {
    std::ofstream file;
    file.open(filename);
    auto increments = lengths_;
    double a0 =  0.529177210903;  // in Angstrom
    for (auto i = 0; i < 3; ++i) {
        increments[i] /= dimensions_[i];
    }

    auto origin = center_ - lengths_/2;
    file << comment << '\n'
         << dimensions_[0] << ' ' << dimensions_[1] << ' ' << dimensions_[2] << "   //nx ny nz\n"
         << ToString::doubleToString(origin[0] * a0, 5, 0, false) << ' '
         << ToString::doubleToString(origin[1] * a0, 5, 0, false) << ' '
         << ToString::doubleToString(origin[2] * a0, 5, 0, false) << "   //Origin of the 3D grid\n"
         << ToString::doubleToString(increments[0]*a0, 6, 0, false) << ' '
         << ToString::doubleToString(increments[1]*a0, 6, 0, false) << ' '
         << ToString::doubleToString(increments[2]*a0, 6, 0, false)
         << "   //x increment, y inc, z inc/ grid(x(y(z)))\n";

    double totalWeightd = 0.0;
    for (auto &weight : data_){
        totalWeightd += double(weight);
    }

    int l = 0;
    for (IndexType i = 0; i < dimensions_[0]; ++i) {
        for (IndexType j = 0; j < dimensions_[1]; ++j) {
            for (IndexType k = 0; k < dimensions_[2]; ++k) {
                file << ToString::doubleToString(
                        double(data_[VoxelCube::index(i, j, k)])
                        / (increments.prod()) / totalWeightd, 11, 0, false);
                if ((l + 1) % 5 == 0) {
                    file << '\n';
                } else {
                    file << ' ';
                }
                ++l;
            }
        }
    }
}


namespace YAML {
    Node convert<VoxelCube>::encode(const VoxelCube &rhs) {
        Node node;
        node["smoothed"] = rhs.smoothQ_;
        node["dimensions"] = rhs.getDimensions();
        node["lengths"] = rhs.getLengths();
        node["center"] = rhs.getCenter();
        node["data"] = rhs.getData();
        return node;
    }

    bool convert<VoxelCube>::decode(const Node &node, VoxelCube &rhs) {
        if (!node.IsMap())
            return false;

        bool smoothed = false;
        if(node["smoothed"])
            smoothed = node["smoothed"].as<bool>();
        auto dimensions = node["dimensions"].as<Eigen::Matrix<VoxelCube::IndexType, 3, 1>>();
        auto lengths = node["lengths"].as<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>();
        auto center = node["center"].as<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>();
        rhs = VoxelCube(dimensions, lengths, center, smoothed);
        rhs.setData(node["data"].as<std::vector<VoxelCube::VolumeDataType>>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const VoxelCube &rhs) {
        out << BeginMap
            << Key << "smoothed" << Value << rhs.smoothQ_
            << Key << "dimensions" << Value << rhs.getDimensions()
            << Key << "lengths" << Value << rhs.getLengths()
            << Key << "center" << Value << rhs.getCenter()
            << Key << "data" << Value << Flow << BeginSeq;
        for (const auto &i : rhs.getData()) {
            out << int(i);
        }
        out << EndSeq << EndMap;
        return out;
    }
}
