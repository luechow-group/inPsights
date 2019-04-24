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