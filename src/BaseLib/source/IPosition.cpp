/* Copyright (C) 2019 Leonard Reuter.
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
