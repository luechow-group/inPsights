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

template <typename T> //uint_8t or uint16_t
class VoxelCube {
    using IndexType = dualmc::QuadIndexType;
    using VertexComponentsType = dualmc::VertexComponentsType;
public:
    explicit VoxelCube(
            IndexType dimension = 128,
            VertexComponentsType length = 8 * ConversionFactors::angstrom2bohr)
            :
            dimension(dimension),
            length(length), halfLength(length/2.0), inverseDimension(1.0/(dimension-1)),
            data(static_cast<unsigned long>(dimension*dimension*dimension))
    {
        assert(dimension > 3 && "The cube dimensions must be greater than 3.");
        assert(length > 0 && "Length must be positive.");
    }

    inline long index(IndexType x, IndexType y, IndexType z){
        return x + dimension * (y + dimension*z);
    }

    void add(const Eigen::Vector3d& pos, IndexType weight = 1) {
        if(pos.minCoeff() >= -halfLength && pos.maxCoeff() < halfLength ) {
            auto x = IndexType((pos[0] + halfLength) / length * dimension);
            auto y = IndexType((pos[1] + halfLength) / length * dimension);
            auto z = IndexType((pos[2] + halfLength) / length * dimension);

            data[index(x,y,z)] += weight;
        }
    };

    void shiftDualMCResults(std::vector<dualmc::Vertex>& vertices){
        for(auto & v : vertices){
            v.x = (v.x * inverseDimension - offset)*length;
            v.y = (v.y * inverseDimension - offset)*length;
            v.z = (v.z * inverseDimension - offset)*length;
        }
    }

    /*void shiftDualMCResults(std::vector<Vertex>& vertices){
        for(auto & v : vertices){
            //v.position /= dimension; // TODO use invDim ??
            v.position *= inverseDimension;
            v.position.array() -= offset;
            v.position *= length;
        }
    }*/

    IndexType dimension;
    VertexComponentsType length, halfLength, inverseDimension;
    std::vector<T> data;
    const VertexComponentsType offset = 0.5;
};

#endif //INPSIGHTS_VOXELCUBE_H
