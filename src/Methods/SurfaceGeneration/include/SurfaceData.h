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
