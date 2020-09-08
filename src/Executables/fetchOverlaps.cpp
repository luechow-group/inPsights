// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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