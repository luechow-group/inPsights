// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
