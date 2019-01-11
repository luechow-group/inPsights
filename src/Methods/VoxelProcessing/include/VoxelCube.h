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

class VoxelCube {
public:
    using IndexType = dualmc::QuadIndexType;
    using VertexComponentsType = dualmc::VertexComponentsType;

    explicit VoxelCube(
            IndexType dimension = 32,
            VertexComponentsType length = VertexComponentsType(8 * ConversionFactors::angstrom2bohr),
            const Eigen::Matrix<VertexComponentsType,3,1>& origin = {0,0,0});

    inline long index(IndexType x, IndexType y, IndexType z);

    void add(const Eigen::Vector3d& pos, IndexType weight = 1);

    void shiftDualMCResults(std::vector<dualmc::Vertex>& vertices);

    IndexType dimension;
    VertexComponentsType length, halfLength, inverseDimension;
    std::vector<uint16_t > data;
    Eigen::Matrix<VertexComponentsType,3,1> origin;
    static constexpr VertexComponentsType offset = 0.5;
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
