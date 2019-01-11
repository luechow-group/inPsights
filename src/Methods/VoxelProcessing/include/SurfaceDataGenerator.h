//
// Created by Michael Heuer on 2019-01-07.
//

#ifndef INPSIGHTS_SURFACEDATAGENERATOR_H
#define INPSIGHTS_SURFACEDATAGENERATOR_H

#include <Eigen/Core>
#include <DualMC.h>
#include <cstdint>
#include "VoxelCube.h"
#include "Conversion.h"
#include "SurfaceData.h"
#include <Logger.h>
#include <algorithm>
#include <utility>


class SurfaceDataGenerator {
public:
    uint16_t getIsoValue(double volumeThreshold = 0.5, unsigned maxSteps = 10, double eps = 1e-2) {
        auto res = std::minmax_element(cube.data.begin(), cube.data.end());
        auto min = *res.first, max = *res.second;
        auto lower = min, mid = uint16_t(0.2*max), upper = max;

        for (unsigned i = 0; i < maxSteps; ++i) {
            unsigned sumIn = 0, sumOut = 0;
            for (auto v : cube.data) {
                if (v > mid)
                    sumIn += v;
                else
                    sumOut += v;
            }
            auto ratio = double(sumIn)/double(sumIn+sumOut);
            if(std::abs(ratio-volumeThreshold) < eps)
                break;
            else if(ratio < volumeThreshold)
                upper = mid;
            else
                lower = mid;
            mid = (upper+lower) / uint16_t(2);
        }
        return mid;
    }

    explicit SurfaceDataGenerator(VoxelCube volume)
            : cube(volume), dualMcVertices_(), dualMcQuads_() {} //TODO use move semantics?

    SurfaceData computeSurfaceData(double volumeThreshold = 0.5) {
        dualmc::DualMC<uint16_t> builder;

        auto isovalue = getIsoValue(volumeThreshold);

        auto quadSoupQ = false;
        auto manifoldQ = true;
        builder.build(&cube.data.front(),
                      cube.dimension, cube.dimension, cube.dimension,
                      isovalue,
                      manifoldQ, quadSoupQ, dualMcVertices_, dualMcQuads_);

        // correct vertices by voxel cube dimension, offset and length
        cube.shiftDualMCResults(dualMcVertices_);

        Logger::console->info("{0} vertices and {1} quads found for an iso value {2}",
                dualMcVertices_.size(), dualMcQuads_.size(), isovalue);


        SurfaceData surfaceData;
        surfaceData.triangles = Conversion::quadsToTriangles(dualMcQuads_);
        surfaceData.vertices = Conversion::convertVertices(dualMcVertices_);
        
        Conversion::calculateVertexNormals(surfaceData.vertices,surfaceData.triangles);

        return surfaceData;
    }

    VoxelCube cube; //TODO template
    std::vector<dualmc::Vertex> dualMcVertices_;
    std::vector<dualmc::Quad> dualMcQuads_;
};

#endif //INPSIGHTS_SURFACEDATAGENERATOR_H
