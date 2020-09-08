// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
