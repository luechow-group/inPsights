//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Group.h>
#include <BestMatchDistanceIdentityClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <BestMatchDistanceDensityBasedClusterer.h>
#include <BestMatchSOAPSimilarityClusterer.h>
#include <MaximaProcessor.h>
#include <GeneralStatistics.h>
#include <algorithm>
#include <utility>
#include <MaximaProcessingSettings.h>
#include <VoxelCubeGeneration.h>
#include <spdlog/spdlog.h>

using namespace YAML;

int main(int argc, char *argv[]) {

    std::string inputFilename = argv[1];
    YAML::Node inputYaml = YAML::LoadFile(argv[1]);
    YAML::Emitter emitter;
    emitter << inputYaml;

    YAML::Emitter outputYaml;
    spdlog::info("Executable: {}", argv[0]);
    spdlog::info("Input file {}:\n{}", inputFilename, emitter.c_str());

    // Apply settings from inputYaml
    Settings::MaximaProcessing maximaProcessingSettings(inputYaml);

    // Read maxima and samples
    Group maxima;
    std::vector<Sample> samples;
    RawDataReader reader(maxima, samples);
    reader.read(maximaProcessingSettings.binaryFileBasename(), maximaProcessingSettings.samplesToAnalyze());
    auto atoms = reader.getAtoms();

    spdlog::info("number of inital refs {}", maxima.size());
    auto results = GeneralStatistics::calculate(maxima, samples, atoms);

    auto valueStandardError = results.valueStats_.standardError()(0,0);

    // write used settings
    YAML::Node usedSettings, usedClusteringSettings;
    maximaProcessingSettings.appendToNode(usedSettings);

    for(auto node : inputYaml["Clustering"]){
        auto methodName = node.first.as<std::string>();
        spdlog::info("Starting {}...", methodName);

        switch (IClusterer::typeFromString(methodName)) {
            case IClusterer::Type::BestMatchDistanceIdentityClusterer: {
                BestMatchDistanceIdentityClusterer::settings =
                        Settings::BestMatchDistanceIdentityClusterer(inputYaml["Clustering"]);

                BestMatchDistanceIdentityClusterer bestMatchDistanceIdentityClusterer(samples);
                if(!inputYaml["BestMatchDistanceIdentityClusterer"]["identityValueIncrement"])
                    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = valueStandardError*1e-4;
                bestMatchDistanceIdentityClusterer.cluster(maxima);

                BestMatchDistanceIdentityClusterer::settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::BestMatchDistanceSimilarityClusterer: {
                BestMatchDistanceSimilarityClusterer::settings =
                        Settings::BestMatchDistanceSimilarityClusterer(inputYaml["Clustering"]);
                BestMatchDistanceSimilarityClusterer bestMatchDistanceSimilarityClusterer(samples);
                if(!inputYaml["BestMatchDistanceSimilarityClusterer"]["similarityValueIncrement"])
                    BestMatchDistanceSimilarityClusterer::settings.similarityValueIncrement = valueStandardError*1e-2;
                bestMatchDistanceSimilarityClusterer.cluster(maxima);

                BestMatchDistanceSimilarityClusterer::settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::BestMatchDistanceDensityBasedClusterer: {
                BestMatchDistanceSimilarityClusterer::settings =
                        Settings::BestMatchDistanceSimilarityClusterer(inputYaml["Clustering"]);
                BestMatchDistanceDensityBasedClusterer bestMatchDistanceDensityBasedClusterer(samples);
                bestMatchDistanceDensityBasedClusterer.cluster(maxima);

                BestMatchDistanceDensityBasedClusterer::settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::BestMatchSOAPSimilarityClusterer: {
                BestMatchSOAPSimilarityClusterer::settings =
                        Settings::BestMatchSOAPSimilarityClusterer(inputYaml["Clustering"]);
                BestMatchSOAPSimilarityClusterer bestMatchSOAPSimilarityClusterer(atoms, samples);
                bestMatchSOAPSimilarityClusterer.cluster(maxima);

                BestMatchSOAPSimilarityClusterer::settings.appendToNode(usedClusteringSettings);
                break;
            }
            default:
                throw std::invalid_argument("Unknown clustering method \"" + methodName +  "\".");
        }
        spdlog::info("number of elements after {}: {}", methodName, maxima.size());
    }
    maxima.sortAll();

    usedSettings["Clustering"] = usedClusteringSettings;

    VoxelCubeGeneration::settings = Settings::VoxelCubeGeneration(inputYaml);
    VoxelCubeGeneration::settings.appendToNode(usedSettings);

    outputYaml << BeginDoc << Comment("used settings") << usedSettings << EndDoc;

    // write results
    outputYaml << BeginDoc << BeginMap
               << Key << "Atoms" << Value << atoms << Comment("[a0]")
               << Key << "NSamples" << Value << samples.size()
               << Key << "OverallResults" << Value << results;
    spdlog::info("Calculating statistics...");

    MaximaProcessor maximaProcessor(outputYaml, samples,atoms);
    maximaProcessor.calculateStatistics(maxima);

    outputYaml << EndMap << EndDoc;

    std::string resultsFilename = maximaProcessingSettings.binaryFileBasename() + ".yml";
    spdlog::info("Writing results into file \"{}\"", resultsFilename);

    std::ofstream yamlFile(resultsFilename);
    yamlFile << outputYaml.c_str();
    yamlFile.close();

    spdlog::info("Done! Bye bye.");

    return 0;

    /*TODO
     *  - has an averagedStructure
     *  - similarity attribute (enum) - spatially|permutationally|rotationally + similar|identical,
     * - use averaged structure for permutation sort, store best match permutation to add energies
     * - make single value statistics class
     * - test naive std
     * - choice of function value increment
     * - validate that ring-like clusters are ordered correctly
     * - use global similarity for permutation sorting
     * - split identity sort into batches that can be compared in parallel using OpenMP
     * - improve spinSpecificHungarian (low priority)
     * */
};
