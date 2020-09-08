// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SURFACEDATA_H
#define INPSIGHTS_SURFACEDATA_H

#include "Triangle.h"
#include "Vertex.h"
#include <vector>

struct SurfaceData{

    SurfaceData() = default;
    SurfaceData(const std::vector<Vertex>& vertices, const std::vector<Triangle>& triangles)
    : vertices(vertices), triangles(triangles){}

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<SurfaceData> {
        static Node encode(const SurfaceData &rhs);
        static bool decode(const Node &node, SurfaceData &rhs);
    };
    Emitter &operator<<(Emitter &out, const SurfaceData &rhs) ;
}


#endif //INPSIGHTS_SURFACEDATA_H
