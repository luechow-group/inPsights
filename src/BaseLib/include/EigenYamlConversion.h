//
// Created by Michael Heuer on 21.06.18.
//

#ifndef INPSIGHTS_EIGENYAMLCONVERSION_H
#define INPSIGHTS_EIGENYAMLCONVERSION_H

#include<Eigen/Core>

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<>
    struct convert<Eigen::Vector3i> {
        static Node encode(const Eigen::Vector3i& rhs);
        static bool decode(const Node& node, Eigen::Vector3i& rhs);
    };
    Emitter& operator<< (Emitter& out, const Eigen::Vector3i& rhs);

    template<>
    struct convert<Eigen::Vector3f> {
        static Node encode(const Eigen::Vector3f& rhs);
        static bool decode(const Node& node, Eigen::Vector3f& rhs);
    };
    Emitter& operator<< (Emitter& out, const Eigen::Vector3f& rhs);

    template<>
    struct convert<Eigen::Vector3d> {
        static Node encode(const Eigen::Vector3d& rhs);
        static bool decode(const Node& node, Eigen::Vector3d& rhs);
    };
    Emitter& operator<< (Emitter& out, const Eigen::Vector3d& rhs);

    template<>
    struct convert<Eigen::VectorXd> {
        static Node encode(const Eigen::VectorXd& rhs);
        static bool decode(const Node& node, Eigen::VectorXd& rhs);
    };
    Emitter& operator<< (Emitter& out, const Eigen::VectorXd& rhs);

}

#endif //INPSIGHTS_EIGENYAMLCONVERSION_H
