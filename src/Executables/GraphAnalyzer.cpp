/* Copyright (C) 2019 Michael Heuer.
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

#include <string>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "Enumerate.h"
#include <iostream>
#include <fstream>
#include "ParticlesVector.h"
#include "MolecularGeometry.h"
#include "BestMatchDistance.h"
#include <iomanip>

struct Cluster {
    long id;
    std::vector<long> yamlIds;
    double weight, minValue;
    std::vector<double> spinCorr;
    std::vector<ElectronsVector> representativeStructures;
    bool foundQ;
};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.yml" << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];

    std::vector<std::vector<Cluster>> clusterData;

    auto results = YAML::LoadAllFromFile(inputFilename)[1];
    auto clustersNodes = results["Clusters"];

    std::vector<ElectronsVector> structures;

    std::size_t limit = clustersNodes.size();
    for (std::size_t k = 0; k < limit; ++k) {
        auto structureNode = clustersNodes[k]["Structures"];
        for (auto structure = structureNode.begin(); structure != structureNode.end(); ++structure) {
            structures.emplace_back((*structure).as<ElectronsVector>());
        }
    }

    double radius = 0.15;

    std::ofstream file("out.txt");

    for (std::size_t i = 0; i < structures.size()-1; ++i) {
        for (std::size_t j = i+1; j < structures.size(); ++j) {
            auto[norm, perm] = BestMatch::Distance::compare(
                    structures[i].positionsVector(),
                    structures[j].positionsVector());
            //if(norm < radius)
                file << "{" << i << "," << j << "," <<  norm << "},";
        }
    }
    file.close();

}