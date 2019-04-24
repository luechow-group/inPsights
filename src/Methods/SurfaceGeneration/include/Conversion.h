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

#ifndef INPSIGHTS_CONVERSION_H
#define INPSIGHTS_CONVERSION_H

#include "Vertex.h"
#include "Triangle.h"
#include <DualMC.h>
#include <vector>

namespace Conversion {
    std::vector<Vertex> convertVertices(const std::vector<dualmc::Vertex>& dualMcVertices);

    std::pair<Triangle, Triangle> quadToTrianglePair(const dualmc::Quad &quad);

    std::vector<Triangle> quadsToTriangles(const std::vector<dualmc::Quad> &quads);

    Eigen::Vector3f calculateFaceNormal(
            const Vertex &v1,
            const Vertex &v2,
            const Vertex &v3);

    void calculateVertexNormals(
            std::vector<Vertex> &vertices,
            const std::vector<Triangle> &triangles);
}

#endif //INPSIGHTS_CONVERSION_H
