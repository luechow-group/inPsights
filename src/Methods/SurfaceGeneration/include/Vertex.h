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

#ifndef INPSIGHTS_VERTEX_H
#define INPSIGHTS_VERTEX_H

#include <Eigen/Core>

struct Vertex {
    Vertex() = default;

    Vertex(Vertex const &v) = default;

    Vertex(const Eigen::Vector3f& position);

    Eigen::Vector3f position, normal;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<Vertex> {
        static Node encode(const Vertex &rhs);
        static bool decode(const Node &node, Vertex &rhs);
    };
    Emitter &operator<<(Emitter &out, const Vertex &rhs) ;
}

#endif //INPSIGHTS_VERTEX_H