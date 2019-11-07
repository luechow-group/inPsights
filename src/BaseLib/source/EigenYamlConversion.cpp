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

#include "EigenYamlConversion.h"
#include <yaml-cpp/yaml.h>

namespace YAML {
    Node convert<Eigen::Vector3i>::encode(const Eigen::Vector3i &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        return node;
    }
    bool convert<Eigen::Vector3i>::decode(const Node &node, Eigen::Vector3i &rhs) {
        if (!node.IsSequence() || node.size() != 3)
            return false;

        rhs[0] = node[0].as<int>();
        rhs[1] = node[1].as<int>();
        rhs[2] = node[2].as<int>();
        return true;
    }
    Emitter &operator<<(Emitter &out, const Eigen::Vector3i &v) {
        out << Flow << BeginSeq << v[0] << v[1] << v[2] << EndSeq;
        return out;
    }

    Node convert<Eigen::Vector3f>::encode(const Eigen::Vector3f &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        return node;
    }
    bool convert<Eigen::Vector3f>::decode(const Node &node, Eigen::Vector3f &rhs) {
        if (!node.IsSequence() || node.size() != 3)
            return false;

        rhs[0] = node[0].as<float>();
        rhs[1] = node[1].as<float>();
        rhs[2] = node[2].as<float>();
        return true;
    }
    Emitter &operator<<(Emitter &out, const Eigen::Vector3f &v) {
        out << Flow << BeginSeq << v[0] << v[1] << v[2] << EndSeq;
        return out;
    }

    Node convert<Eigen::Vector3d>::encode(const Eigen::Vector3d &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        return node;
    }
    bool convert<Eigen::Vector3d>::decode(const Node &node, Eigen::Vector3d &rhs) {
        if (!node.IsSequence() || node.size() != 3)
            return false;

        rhs[0] = node[0].as<double>();
        rhs[1] = node[1].as<double>();
        rhs[2] = node[2].as<double>();
        return true;
    }
    Emitter &operator<<(Emitter &out, const Eigen::Vector3d &v) {
        out << Flow << BeginSeq << v[0] << v[1] << v[2] << EndSeq;
        return out;
    }

    Node convert<Eigen::VectorXd>::encode(const Eigen::VectorXd &rhs) {
        Node node;

        for (Eigen::Index i = 0; i < rhs.size(); ++i)
            node.push_back(rhs[i]);

        return node;
    }
    bool convert<Eigen::VectorXd>::decode(const Node &node, Eigen::VectorXd &rhs) {
        if (!node.IsSequence())
            return false;

        rhs.resize(node.size());

        for (size_t i = 0; i < node.size(); ++i)
            rhs[i] = node[i].as<double>();

        return true;
    }
    Emitter &operator<<(Emitter &out, const Eigen::VectorXd &v) {
        out << Flow << BeginSeq;
        for (Eigen::Index i = 0; i < v.size(); ++i) out << v[i];
        out << EndSeq;
        return out;
    }

    Node convert<Eigen::MatrixXd>::encode(const Eigen::MatrixXd &rhs) {
        Node node;

        for (Eigen::Index i = 0; i < rhs.rows(); ++i) {
            node["."] = 0; // Trick to get a map
            node.remove(".");
            for (Eigen::Index j = 0; j < rhs.cols(); ++j) {
                node[i][j] = rhs(i,j);
                node[i]["."] = 0; // Trick to get a map
                node[i].remove(".");
            }
        }

        return node;
    }
    bool convert<Eigen::MatrixXd>::decode(const Node &node, Eigen::MatrixXd &rhs) {

        if (!node.IsMap() || !node[0].IsMap())
            return false;

        rhs.resize(node.size(),node[0].size());

        for (Eigen::Index i = 0; i < rhs.rows(); ++i){
            for (Eigen::Index j = 0; j < rhs.cols(); ++j){
                rhs(i,j) = node[i][j].as<double>();
            }
        }

        return true;
    }
    Emitter &operator<<(Emitter &out, const Eigen::MatrixXd &rhs) {

        out << BeginMap;
        for (Eigen::Index i = 0; i < rhs.rows(); ++i) {
            out << Key << i << Value << BeginMap;
            for (Eigen::Index j = 0; j < rhs.cols(); ++j) {
                out << Key << j << Value << rhs(i,j);
            }
            out << EndMap;
        }
        out << EndMap;
        return out;
    }

    Node convert<TriangularMatrixXd>::encode(
            const TriangularMatrixXd &rhs) {
        Node node;

        Eigen::MatrixXd mat(rhs);

        for (Eigen::Index i = 0; i < mat.rows()-1; ++i) {
            node["."] = 0; // Trick to get a map
            node.remove(".");
            for (Eigen::Index j = i+1; j < mat.cols(); ++j) {
                node[i][j] = mat(i,j);
                node[i]["."] = 0; // Trick to get a map
                node[i].remove(".");
            }
        }

        return node;
    }
    bool convert<TriangularMatrixXd>::decode(const Node &node, TriangularMatrixXd &rhs) {

        if (!node.IsMap() || !node[0].IsMap())
            return false;

        auto& mat = rhs.nestedExpression();
        mat = Eigen::MatrixXd::Zero(node.size()+1,node[0].size()+1);

        for (Eigen::Index i = 0; i < mat.rows()-1; ++i){
            for (Eigen::Index j = i+1; j < mat.cols(); ++j){
                mat(i,j) = node[i][j].as<double>();
            }
        }

        return true;
    }
    Emitter &operator<<(Emitter &out, const TriangularMatrixXd &rhs) {

        auto& mat = rhs.nestedExpression();

        out << BeginMap;
        for (Eigen::Index i = 0; i < mat.rows()-1; ++i) {
            out << Key << i << Value << BeginMap;
            for (Eigen::Index j = i+1; j < mat.cols(); ++j) {
                out << Flow << Key << j << Value << mat(i,j);
            }
            out << EndMap;
        }
        out << EndMap;

        return out;
    }
}