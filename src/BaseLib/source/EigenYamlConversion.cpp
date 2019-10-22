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

        for (long i = 0; i < rhs.size(); ++i)
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
        for (long i = 0; i < v.size(); ++i) out << v[i];
        out << EndSeq;
        return out;
    }
}