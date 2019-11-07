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

using TriangularMatrixXd = Eigen::MatrixXd::SelfAdjointViewReturnType<Eigen::Upper>::Type;//Eigen::MatrixBase<Eigen::MatrixXd>::SelfAdjointViewReturnType<Eigen::Upper>::Type;

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

    template<>
    struct convert<Eigen::MatrixXd> {
        static Node encode(const Eigen::MatrixXd& rhs);
        static bool decode(const Node& node, Eigen::MatrixXd& rhs);
    };
    Emitter& operator<< (Emitter& out, const Eigen::MatrixXd& rhs);

    template<>
    struct convert<TriangularMatrixXd> {
        static Node encode(const TriangularMatrixXd& rhs);
        static bool decode(const Node& node, TriangularMatrixXd& rhs);
    };
    Emitter& operator<< (Emitter& out, const TriangularMatrixXd& rhs);
}

#endif //INPSIGHTS_EIGENYAMLCONVERSION_H
