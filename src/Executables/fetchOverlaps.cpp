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

#include <iostream>
#include <yaml-cpp/yaml.h>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.yml " << "clusterID" << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];
    int clusterId = std::atoi(argv[2]);

    YAML::Node doc = YAML::LoadAllFromFile(argv[1])[1];

    std::cout << doc["Clusters"][clusterId]["SedOverlaps"] << std::endl;
}