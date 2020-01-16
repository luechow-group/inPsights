/* Copyright 2020 Michael Heuer.
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
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <ParticlesVector.h>
#include <SpinCorrelationValueDistribution.h>
#include <fstream>

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " input.yml " << " oneSidedBinCount" << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];

    int oneSidedBinCount = 12;
    if (argc == 3)
        oneSidedBinCount = std::atoi(argv[2]);
    std::cout << oneSidedBinCount << std::endl;

    YAML::Node doc = YAML::LoadFile(inputFilename);
    const YAML::Node& filenamesNode = doc["files"];

    YAML::Emitter outputYamlDoc;


    Eigen::VectorXd globalStats = Eigen::VectorXd::Zero(2*oneSidedBinCount+1);

    if(filenamesNode.IsSequence()) {

        outputYamlDoc << YAML::BeginMap;
        for (YAML::const_iterator it = filenamesNode.begin(); it != filenamesNode.end(); ++it) {
            auto filename = (*it).as<std::string>();
            std::cout << filename << std::endl;

            YAML::Node molDoc = YAML::LoadAllFromFile(filename)[1];

            SpinCorrelationValueDistribution distribution(oneSidedBinCount);
            Eigen::VectorXd molecularStats = Eigen::VectorXd::Zero(2*oneSidedBinCount+1);

            //auto nElectrons = molDoc["Clusters"][0]["Structures"][0].as<ElectronsVector>().numberOfEntities();

            auto nClusters = static_cast<int>(molDoc["Clusters"].size());
            for (int clusterId = 0; clusterId < nClusters; ++clusterId) {
                auto spinCorrStats = molDoc["Clusters"][clusterId]["SpinCorrelations"].as<TriangularMatrixStatistics>();
                distribution.addSpinStatistic(spinCorrStats);
            }
            molecularStats += distribution.getHistogramVector();

            outputYamlDoc << YAML::Key << filename << YAML::Value << molecularStats;
            std::cout << molecularStats.transpose()  << std::endl;

            globalStats += molecularStats;
        }

    } else {
        spdlog::error("The input file key \"files\" must contain a sequence.");
    }

    Eigen::VectorXd globalStatsNormalized= globalStats/globalStats.sum();
    outputYamlDoc
    << YAML::Key << "Global" << YAML::Value << globalStats
    << YAML::Key << "GlobalNormalized" << YAML::Value << globalStatsNormalized << YAML::EndMap;

    std::cout << std::endl << globalStats.transpose() << std::endl;
    std::cout << globalStatsNormalized.transpose() << std::endl;

    // write output yaml
    std::string resultsBaseName = inputFilename.substr(0,inputFilename.find('.')) + "-out";
    std::string resultsFilename = resultsBaseName + ".yml";

    if (std::ifstream(resultsFilename).good()){
        for (int i = 1; i < 100; i++){
            resultsFilename = resultsBaseName + "-" + std::to_string(i) + ".yml";
            if (not std::ifstream(resultsFilename).good()){
                break;
            }
        }
    }

    std::ofstream yamlFile(resultsFilename);
    yamlFile << outputYamlDoc.c_str();
    yamlFile.close();
    std::cout << "output wrote to file " << resultsFilename << std::endl;
}