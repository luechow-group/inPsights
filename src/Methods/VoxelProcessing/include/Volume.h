//
// Created by Michael Heuer on 2019-01-08.
//

#ifndef INPSIGHTS_VOLUME_H
#define INPSIGHTS_VOLUME_H

#include <cstdint>
#include <vector>
#include <NaturalConstants.h>
#include <Eigen/Core>

template <typename T> // uint16_t or uint_6t
struct Volume {
    Volume(double length = 8*ConversionFactors::angstrom2bohr, uint8_t numberOfVoxels = 32)
            : length(length), halfLength(length/2.0),
            dimX(numberOfVoxels), dimY(numberOfVoxels), dimZ(numberOfVoxels),
            invDimX(1.0/dimX), invDimY(1.0/dimY), invDimZ(1.0/dimZ),
            data(numberOfVoxels*numberOfVoxels*numberOfVoxels)
    {}

    double length, halfLength;
    uint8_t dimX, dimY, dimZ;

    double invDimX, invDimY, invDimZ;
    std::vector<T> data;

    inline long index(uint8_t x, uint8_t y, uint8_t z){
        return x*dimY*dimZ + y*dimZ + z;
    }

    void add(const Eigen::Vector3d& pos, uint8_t weight = 1) {
        if(pos.minCoeff() >= -halfLength && pos.maxCoeff() < halfLength ) {
            uint8_t x = uint8_t((pos[0] + halfLength) / length * invDimX);
            uint8_t y = uint8_t((pos[1] + halfLength) / length * invDimY);
            uint8_t z = uint8_t((pos[2] + halfLength) / length * invDimZ);

            data[index(x,y,z)] += weight;
        }
    };
};

#endif //INPSIGHTS_VOLUME_H
