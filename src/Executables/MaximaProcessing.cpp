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
#include <SOAPSettings.h>
#include <AsciiArt.h>

using namespace YAML;

void validateClusteringSettings(const YAML::Node &inputYaml) {

    auto clusteringNode = inputYaml["Clustering"];

    for (auto node : clusteringNode) {
        auto methodName = node.first.as<std::string>();

        switch (IClusterer::typeFromString(methodName)) {
            case IClusterer::Type::BestMatchDistanceIdentityClusterer: {
                BestMatchDistanceIdentityClusterer::settings
                        = Settings::BestMatchDistanceIdentityClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::BestMatchDistanceSimilarityClusterer: {
                BestMatchDistanceSimilarityClusterer::settings
                        = Settings::BestMatchDistanceSimilarityClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::BestMatchDistanceDensityBasedClusterer: {
                BestMatchDistanceDensityBasedClusterer::settings
                        = Settings::BestMatchDistanceDensityBasedClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::BestMatchSOAPSimilarityClusterer: {
                BestMatchSOAPSimilarityClusterer::settings
                        = Settings::BestMatchSOAPSimilarityClusterer(clusteringNode);

                auto soapSettings = node.second["SOAP"];
                if (soapSettings[SOAP::General::settings.name()])
                    SOAP::General::settings = Settings::SOAP::General(soapSettings);
                if (soapSettings[SOAP::Radial::settings.name()])
                    SOAP::Radial::settings = Settings::SOAP::Radial(soapSettings);
                if (soapSettings[SOAP::Angular::settings.name()])
                    SOAP::Angular::settings = Settings::SOAP::Angular(soapSettings);
                if (soapSettings[SOAP::Cutoff::settings.name()])
                    SOAP::Cutoff::settings = Settings::SOAP::Cutoff(soapSettings);
                break;
            }
            default:
                throw std::invalid_argument("Unknown clustering method \"" + methodName + "\".");
        }
    }
}

void validateInput(const YAML::Node &inputYaml) {
    spdlog::set_level(spdlog::level::off);

    Settings::MaximaProcessing maximaProcessingSettings(inputYaml);
    validateClusteringSettings(inputYaml);
    VoxelCubeGeneration::settings = Settings::VoxelCubeGeneration(inputYaml);

    spdlog::set_level(spdlog::level::info);
    spdlog::info("Input is valid.");
}

int main(int argc, char *argv[]) {

    std::string inputFilename = argv[1];
    YAML::Node inputYaml = YAML::LoadFile(argv[1]);
    YAML::Emitter emitter;
    emitter << inputYaml;

    YAML::Emitter outputYaml;
    spdlog::info("Executable: {}", argv[0]);
    spdlog::info("Validating input from file {}:\n{}...", inputFilename, emitter.c_str());
    validateInput(inputYaml);


    // Apply settings from inputYaml
    MaximaProcessing::settings = Settings::MaximaProcessing(inputYaml);

    // Read maxima and samples
    Group maxima;
    std::vector<Sample> samples;
    RawDataReader reader(maxima, samples);
    reader.read(MaximaProcessing::settings.binaryFileBasename(), MaximaProcessing::settings.samplesToAnalyze());
    auto atoms = reader.getAtoms();

    spdlog::info("number of inital refs {}", maxima.size());

    // calculate general statistics
    auto results = GeneralStatistics::calculate(maxima, samples, atoms);
    auto valueStandardError = results.valueStats_.standardError()(0, 0);

    // write used settings
    YAML::Node usedSettings, usedClusteringSettings;
    MaximaProcessing::settings.appendToNode(usedSettings);

    auto clusteringNode = inputYaml["Clustering"];

    for (auto node : clusteringNode) {
        auto methodName = node.first.as<std::string>();
        spdlog::info("Starting {}...", methodName);

        if (usedClusteringSettings[methodName])
            spdlog::warn("Method \"{}\" is being applied multiple times! Overwriting old settings...", methodName);

        switch (IClusterer::typeFromString(methodName)) {
            case IClusterer::Type::BestMatchDistanceIdentityClusterer: {
                auto &settings = BestMatchDistanceIdentityClusterer::settings;

                settings = Settings::BestMatchDistanceIdentityClusterer(clusteringNode);

                BestMatchDistanceIdentityClusterer bestMatchDistanceIdentityClusterer(samples);
                if (!clusteringNode[settings.name()][settings.identityValueIncrement.name()])
                    settings.identityValueIncrement = valueStandardError * 1e-4;

                bestMatchDistanceIdentityClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::BestMatchDistanceSimilarityClusterer: {
                auto &settings = BestMatchDistanceSimilarityClusterer::settings;

                settings = Settings::BestMatchDistanceSimilarityClusterer(clusteringNode);

                BestMatchDistanceSimilarityClusterer bestMatchDistanceSimilarityClusterer(samples);
                if (!clusteringNode[settings.name()][settings.similarityValueIncrement.name()])
                    settings.similarityValueIncrement = valueStandardError * 1e-2;
                bestMatchDistanceSimilarityClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::BestMatchDistanceDensityBasedClusterer: {
                auto &settings = BestMatchDistanceDensityBasedClusterer::settings;

                settings = Settings::BestMatchDistanceDensityBasedClusterer(clusteringNode);

                BestMatchDistanceDensityBasedClusterer bestMatchDistanceDensityBasedClusterer(samples);
                bestMatchDistanceDensityBasedClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::BestMatchSOAPSimilarityClusterer: {
                auto &settings = BestMatchSOAPSimilarityClusterer::settings;

                settings = Settings::BestMatchSOAPSimilarityClusterer(clusteringNode);

                auto soapSettings = node.second["SOAP"];
                if (soapSettings[SOAP::General::settings.name()])
                    SOAP::General::settings = Settings::SOAP::General(soapSettings);
                if (soapSettings[SOAP::Radial::settings.name()])
                    SOAP::Radial::settings = Settings::SOAP::Radial(soapSettings);
                if (soapSettings[SOAP::Angular::settings.name()])
                    SOAP::Angular::settings = Settings::SOAP::Angular(soapSettings);
                if (soapSettings[SOAP::Cutoff::settings.name()])
                    SOAP::Cutoff::settings = Settings::SOAP::Cutoff(soapSettings);

                BestMatchSOAPSimilarityClusterer bestMatchSOAPSimilarityClusterer(atoms, samples);
                bestMatchSOAPSimilarityClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                auto usedSoapSettings = usedClusteringSettings[settings.name()]["SOAP"];
                SOAP::General::settings.appendToNode(usedSoapSettings);
                SOAP::Radial::settings.appendToNode(usedSoapSettings);
                SOAP::Angular::settings.appendToNode(usedSoapSettings);
                SOAP::Cutoff::settings.appendToNode(usedSoapSettings);

                break;
            }
            default:
                throw std::invalid_argument("Unknown clustering method \"" + methodName + "\".");
        }
        spdlog::info("number of elements after {}: {}", methodName, maxima.size());
    }
    maxima.sortAll();

    usedSettings["Clustering"] = usedClusteringSettings;

    VoxelCubeGeneration::settings = Settings::VoxelCubeGeneration(inputYaml);
    VoxelCubeGeneration::settings.appendToNode(usedSettings);


    outputYaml << BeginDoc
    << Comment("input from \"" + inputFilename + "\"") << inputYaml << EndDoc;
    outputYaml << BeginDoc
               << BeginMap
               << Comment(AsciiArt::inPsightsLogo)
               << Key << "Atoms" << Value << atoms << Comment("[a0]")
               << Key << "NSamples" << Value << samples.size()
               << Key << "OverallResults" << Value << results;
    spdlog::info("Calculating statistics...");

    MaximaProcessor maximaProcessor(outputYaml, samples, atoms);
    maximaProcessor.calculateStatistics(maxima);

    outputYaml << EndMap << EndDoc;
    outputYaml << BeginDoc << Comment("final settings") << usedSettings << EndDoc;

    std::string resultsFilename = MaximaProcessing::settings.binaryFileBasename() + ".yml";
    spdlog::info("Writing results into file \"{}\"", resultsFilename);

    std::ofstream yamlFile(resultsFilename);
    yamlFile << outputYaml.c_str();
    yamlFile.close();

    spdlog::info("Done! Bye bye.");

    return 0;

    /*TODO
     * - add similarity attribute (enum - spatially|permutationally|rotationally + similar|identical ) to group
     * - make single value statistics class
     * - test naive standard deviatino
     * - test choice of function value increment
     * - validate that ring-like clusters are ordered correctly
     * - split identity sort into batches that can be compared in parallel using OpenMP
     * - improve spinSpecificHungarian (low priority)
     * */
};
