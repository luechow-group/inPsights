// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_TRIANGLE_H
#define INPSIGHTS_TRIANGLE_H

#include <Eigen/Core>

 struct Triangle {
     Triangle() = default;

     Triangle(const Eigen::Vector3i& indices);

     Eigen::Vector3i indices;
 };

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<Triangle> {
        static Node encode(const Triangle &rhs);
        static bool decode(const Node &node, Triangle &rhs);
    };
    Emitter &operator<<(Emitter &out, const Triangle &rhs) ;
}

#endif //INPSIGHTS_TRIANGLE_H
