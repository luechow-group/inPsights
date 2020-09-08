// Copyright (C) 2019 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <IPosition.h>

namespace YAML{
    Eigen::Vector3d decodePosition(const YAML::Node &node, const AtomsVector &nuclei){
        if (node.size() != 1) throw std::invalid_argument("Nodes in positions must have size 1.");

        if (node["atAtom"].IsDefined()){
            return nuclei[node["atAtom"].as<int>()].position();
        }
        else if (node["atCoordinates"].IsDefined()) {
            return node["atCoordinates"].as<Eigen::Vector3d>();
        }
        else if (node["inBetween"].IsDefined()) {
            Eigen::Vector3d position({0.0,0.0,0.0});
            int positionsNumber = 0;
            for (const auto &positionNode : node["inBetween"]){
                position += decodePosition(positionNode, nuclei);
                positionsNumber++;
            }
            return position / positionsNumber;
        }
        else{
            throw std::invalid_argument("Invalid key in positions.");
        }
    }
}
