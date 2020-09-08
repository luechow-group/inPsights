// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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