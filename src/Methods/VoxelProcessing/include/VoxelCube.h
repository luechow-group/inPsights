//
// Created by Michael Heuer on 2019-01-08.
//

#ifndef INPSIGHTS_VOXELCUBE_H
#define INPSIGHTS_VOXELCUBE_H

#include <cstdint>
#include <vector>
#include <NaturalConstants.h>
#include <Eigen/Core>
#include <DualMC.h>
#include <Vertex.h>

template <typename BinType> //uint_8t or uint16_t
class VoxelCube {
public:
    using IndexType = dualmc::QuadIndexType;
    using VertexComponentsType = dualmc::VertexComponentsType;

    explicit VoxelCube(
            IndexType dimension = 32,
            VertexComponentsType length = 8 * ConversionFactors::angstrom2bohr, const Eigen::Matrix<VertexComponentsType,3,1>& origin = {0,0,0})
            :
            dimension(dimension),
            length(length), halfLength(length/2.0), inverseDimension(1.0/(dimension-1)),
            data(static_cast<unsigned long>(dimension*dimension*dimension)),
            origin(origin)
    {
        assert(dimension > 3 && "The cube dimensions must be greater than 3.");
        assert(length > 0 && "Length must be positive.");
    }

    inline long index(IndexType x, IndexType y, IndexType z){
        return x + dimension * (y + dimension*z);
    }

    void add(const Eigen::Vector3d& pos, IndexType weight = 1) {
        if(pos.minCoeff() >= -halfLength && pos.maxCoeff() < halfLength ) {
            auto x = IndexType((pos[0] + halfLength - origin[0]) / length * dimension);
            auto y = IndexType((pos[1] + halfLength - origin[1]) / length * dimension);
            auto z = IndexType((pos[2] + halfLength - origin[2]) / length * dimension);

            data[index(x,y,z)] += weight;
        }
    };

    void shiftDualMCResults(std::vector<dualmc::Vertex>& vertices){
        for(auto & v : vertices){
            v.x = (v.x * inverseDimension - offset + origin[0])*length;
            v.y = (v.y * inverseDimension - offset + origin[1])*length;
            v.z = (v.z * inverseDimension - offset + origin[2])*length;
        }
    }

    IndexType dimension;
    VertexComponentsType length, halfLength, inverseDimension;
    std::vector<BinType> data;
    static constexpr VertexComponentsType offset = 0.5;
    Eigen::Vector3f origin;
};

#include <yaml-cpp/yaml.h>
namespace YAML {
    template<typename BinType> struct convert<VoxelCube<BinType>> {
        static Node encode(const VoxelCube<BinType> &rhs) {
            Node node;
            node["dimension"] = rhs.dimension;
            node["length"] = rhs.length;
            node["data"] = rhs.data;
            return node;
        }

        static bool decode(const Node &node, VoxelCube<BinType> &rhs) {
            if (!node.IsMap())
                return false;

            auto dimension = node["dimension"].as<typename VoxelCube<BinType>::IndexType>();
            auto length = node["length"].as<typename VoxelCube<BinType>::VertexComponentsType>();

            rhs = VoxelCube<BinType>(dimension, length);
            rhs.data = node["data"].as<std::vector<BinType>>();
            return true;
        }
    };

    template<typename BinType>
    Emitter& operator<< (Emitter& out, const VoxelCube<BinType>& rhs){
        out << BeginMap
        << Key << "dimension" << Value << rhs.dimension
        << Key << "length" << Value << rhs.length
        << Key << "data" << Value << Flow << BeginSeq;
        for(const auto& i : rhs.data) {
            out << int(i);
        }
        out << EndSeq << EndMap;
        return out;
    };
}

#endif //INPSIGHTS_VOXELCUBE_H
