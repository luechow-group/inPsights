// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
