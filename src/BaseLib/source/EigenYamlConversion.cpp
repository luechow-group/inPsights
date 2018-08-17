//
// Created by Michael Heuer on 21.06.18.
//

#include "EigenYamlConversion.h"
#include <yaml-cpp/yaml.h>
#include <ToString.h>

namespace YAML {
    Node convert<Eigen::Vector3d>::encode(const Eigen::Vector3d &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        return node;
    }

    bool convert<Eigen::Vector3d>::decode(const Node &node, Eigen::Vector3d &rhs) {
        if (!node.IsSequence() || node.size() != 3) {
            return false;
        }
        rhs[0] = node[0].as<double>();
        rhs[1] = node[1].as<double>();
        rhs[2] = node[2].as<double>();
        return true;
    }

    Emitter &operator<<(Emitter &out, const Eigen::Vector3d &v) {
        out << Flow << BeginSeq << v[0] << v[1] << v[2] << EndSeq;
        return out;
    }
}