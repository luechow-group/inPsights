// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Triangle.h"
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

Triangle::Triangle(const Eigen::Vector3i& indices)
: indices(indices) {}

namespace YAML {
    Node convert<Triangle>::encode(const Triangle &rhs) {
        Node node;
        node = rhs.indices;
        return node;
    }
    
    bool convert<Triangle>::decode(const Node &node, Triangle &rhs) {
        if (!node.IsSequence())
            return false;

        rhs.indices = node.as<Eigen::Vector3i>();
        return true;
    }

    Emitter &operator<<(Emitter &out, const Triangle &rhs) {
        out << rhs.indices;
        return out;
    }
}