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

#include <Vertex.h>
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

 Vertex::Vertex(const Eigen::Vector3f& position)
 : position(position),normal(Eigen::Vector3f::Zero()) {}

namespace YAML {
    Node convert<Vertex>::encode(const Vertex &rhs) {
        Node node;
        node.push_back(rhs.position);
        node.push_back(rhs.normal);
        return node;
    }
    bool convert<Vertex>::decode(const Node &node, Vertex &rhs) {
     if (!node.IsSequence())
         return false;

     rhs.position = node[0].as<Eigen::Vector3f>();
     rhs.normal = node[1].as<Eigen::Vector3f>();
     return true;
    }

    Emitter &operator<<(Emitter &out, const Vertex &rhs) {
        out << Flow << BeginSeq << rhs.position << rhs.normal << EndSeq;
        return out;
    }
}
