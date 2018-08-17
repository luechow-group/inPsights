//
// Created by Michael Heuer on 21.06.18.
//

#ifndef AMOLQCPP_EIGENYAMLCONVERSION_H
#define AMOLQCPP_EIGENYAMLCONVERSION_H

#include<Eigen/Core>

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<>
    struct convert<Eigen::Vector3d> {
        static Node encode(const Eigen::Vector3d& rhs);
        static bool decode(const Node& node, Eigen::Vector3d& rhs);
    };
    Emitter& operator<< (Emitter& out, const Eigen::Vector3d& p);
}

#endif //AMOLQCPP_EIGENYAMLCONVERSION_H
