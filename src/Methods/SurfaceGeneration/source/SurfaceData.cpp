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

#include "SurfaceData.h"
#include <yaml-cpp/yaml.h>

namespace YAML {
    Node convert<SurfaceData>::encode(const SurfaceData &rhs) {
        Node node;

        node["Vertices"] = rhs.vertices;
        node["Triangles"] = rhs.triangles;

        return node;
    }
    bool convert<SurfaceData>::decode(const Node &node, SurfaceData &rhs) {
        if (!node.IsMap()) {
            return false;
        }
        rhs.vertices = node["Vertices"].as<std::vector<Vertex>>();
        rhs.triangles = node["Triangles"].as<std::vector<Triangle>>();

        return true;
    }

    Emitter &operator<<(Emitter &out, const SurfaceData &rhs) {
        out << BeginMap
        << Key << "Vertices" << Value << rhs.vertices
        << Key << "Triangles" << Value << rhs.triangles
        << EndMap;
        return out;
    }
}
