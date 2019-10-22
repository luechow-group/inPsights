/* Copyright (C) 2018-2019 Michael Heuer.
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

#ifndef INPSIGHTS_EIGENYAMLCONVERSION_H
#define INPSIGHTS_EIGENYAMLCONVERSION_H

#include<Eigen/Core>

/*
namespace YAML {
    template<typename Scalar>
    struct convert<Eigen::Matrix<Scalar,3,1>> {
    static Node encode(const Eigen::Matrix<Scalar,3,1> & rhs){
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        return node;
    }

    static bool decode(const Node& nodes, Eigen::Matrix<Scalar,3,1> & rhs){
        if(!nodes.IsSequence())
            return false;

        Eigen::Matrix<Scalar,3,1> vec;
        for (const auto &i : nodes)
            rhs.append(i.as<Scalar>());

        rhs = vec;
        return true;
    }
};

template<typename Scalar>
Emitter& operator<< (Emitter& out, const Eigen::Matrix<Scalar,3,1>& rhs){
    out << YAML::Flow << YAML::BeginSeq;
    out << rhs[0] << rhs[1] << rhs[3];
    out << YAML::EndSeq;
    return out;
};
}*/

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
