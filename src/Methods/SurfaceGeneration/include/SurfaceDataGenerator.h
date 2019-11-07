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

#ifndef INPSIGHTS_SURFACEDATAGENERATOR_H
#define INPSIGHTS_SURFACEDATAGENERATOR_H

#include <DualMC.h>
#include "VoxelCube.h"
#include "SurfaceData.h"

class SurfaceDataGenerator {
public:

    explicit SurfaceDataGenerator(VoxelCube volume);

    SurfaceData computeSurfaceData(double volumeThreshold = 0.5);

    VoxelCube::VolumeDataType getIsoValue(double volumeThreshold = 0.5, unsigned maxSteps = 10, double eps = 1e-2);

    VoxelCube cube_;
    std::vector<dualmc::Vertex> dualMcVertices_;
    std::vector<dualmc::Quad> dualMcQuads_;
};

#endif //INPSIGHTS_SURFACEDATAGENERATOR_H
