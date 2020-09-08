// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <SurfaceDataGenerator.h>
#include "Conversion.h"
#include <spdlog/spdlog.h>

SurfaceDataGenerator::SurfaceDataGenerator(VoxelCube volume)
        : cube_(std::move(volume)), dualMcVertices_(), dualMcQuads_() {} //TODO use move semantics?

SurfaceData SurfaceDataGenerator::computeSurfaceData(double volumeThreshold) {
    dualmc::DualMC<VoxelCube::VolumeDataType> builder;

    auto isovalue = getIsoValue(volumeThreshold);

    auto quadSoupQ = false;
    auto manifoldQ = true;
    builder.build(&cube_.getData().front(),
                  cube_.getDimension(), cube_.getDimension(), cube_.getDimension(),
                  isovalue,
                  manifoldQ, quadSoupQ, dualMcVertices_, dualMcQuads_);

    // correct vertices by voxel cube dimension, offset and length
    cube_.shiftDualMCResults(dualMcVertices_);

    //spdlog::info("{0} vertices and {1} quads found for an iso value {2}",dualMcVertices_.size(), dualMcQuads_.size(), isovalue);

    SurfaceData surfaceData;
    surfaceData.triangles = Conversion::quadsToTriangles(dualMcQuads_);
    surfaceData.vertices = Conversion::convertVertices(dualMcVertices_);

    Conversion::calculateVertexNormals(surfaceData.vertices, surfaceData.triangles);

    return surfaceData;
}


VoxelCube::VolumeDataType SurfaceDataGenerator::getIsoValue(double volumeThreshold, unsigned maxSteps, double eps) {
    auto res = std::minmax_element(cube_.getData().begin(), cube_.getData().end());
    auto min = *res.first, max = *res.second;
    auto lower = min, mid = VoxelCube::VolumeDataType(0.2 * max), upper = max;

    for (unsigned i = 0; i < maxSteps; ++i) {
        unsigned sumIn = 0, sumOut = 0;
        for (auto v : cube_.getData()) {
            if (v > mid)
                sumIn += v;
            else
                sumOut += v;
        }
        auto ratio = double(sumIn) / double(sumIn + sumOut);
        if (std::abs(ratio - volumeThreshold) < eps)
            break;
        else if (ratio < volumeThreshold)
            upper = mid;
        else
            lower = mid;
        mid = (upper + lower) / VoxelCube::VolumeDataType(2);
    }
    return mid;
}
