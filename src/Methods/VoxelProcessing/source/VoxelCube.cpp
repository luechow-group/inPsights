//
// Created by Michael Heuer on 2019-01-11.
//

#include <VoxelCube.h>
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

VoxelCube::VoxelCube(
        IndexType dimension,
        VertexComponentsType length,
        const Eigen::Matrix<VertexComponentsType, 3, 1> &origin)
        :
        dimension(dimension),
        length(length), halfLength(length / 2.0), inverseDimension(1.0 / (dimension - 1)),
        data(static_cast<unsigned long>(dimension * dimension * dimension)),
        origin(origin) {
    assert(dimension > 3 && "The cube dimensions must be greater than 3.");
    assert(length > 0 && "Length must be positive.");
}

inline long VoxelCube::index(IndexType x, IndexType y, IndexType z) {
    return x + dimension * (y + dimension * z);
}

void VoxelCube::add(const Eigen::Vector3d &pos, IndexType weight) {
    if (pos.minCoeff() >= -halfLength && pos.maxCoeff() < halfLength) {
        auto x = IndexType((pos[0] + halfLength - origin[0]) / length * dimension);
        auto y = IndexType((pos[1] + halfLength - origin[1]) / length * dimension);
        auto z = IndexType((pos[2] + halfLength - origin[2]) / length * dimension);

        data[index(x, y, z)] += weight;
    }
};

void VoxelCube::shiftDualMCResults(std::vector<dualmc::Vertex> &vertices) {
    for (auto &v : vertices) {
        v.x = (v.x * inverseDimension - offset) * length + origin[0];
        v.y = (v.y * inverseDimension - offset) * length + origin[1];
        v.z = (v.z * inverseDimension - offset) * length + origin[2];
    }
}

namespace YAML {
    Node convert<VoxelCube>::encode(const VoxelCube &rhs) {
        Node node;
        node["dimension"] = rhs.dimension;
        node["length"] = rhs.length;
        node["origin"] = rhs.origin;
        node["data"] = rhs.data;
        return node;
    }

    bool convert<VoxelCube>::decode(const Node &node, VoxelCube &rhs) {
        if (!node.IsMap())
            return false;

        auto dimension = node["dimension"].as<VoxelCube::IndexType>();
        auto length = node["length"].as<VoxelCube::VertexComponentsType>();
        auto origin = node["origin"].as<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>();
        rhs = VoxelCube(dimension, length, origin);
        rhs.data = node["data"].as<std::vector<uint16_t >>();
        return true;
    }

    Emitter &operator<<(Emitter &out, const VoxelCube &rhs) {
        out << BeginMap
            << Key << "dimension" << Value << rhs.dimension
            << Key << "length" << Value << rhs.length
            << Key << "origin" << Value << rhs.origin
            << Key << "data" << Value << Flow << BeginSeq;
        for (const auto &i : rhs.data) {
            out << int(i);
        }
        out << EndSeq << EndMap;
        return out;
    }
}
