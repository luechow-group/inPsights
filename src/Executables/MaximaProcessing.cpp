//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Group.h>
#include <IdentityClusterer.h>
#include <DistanceClusterer.h>
#include <DensityBasedClusterer.h>
#include <SOAPClusterer.h>
#include <ReferencePositionsClusterer.h>
#include <MaximaProcessor.h>
#include <GeneralStatistics.h>
#include <algorithm>
#include <utility>
#include <MaximaProcessingSettings.h>
#include <VoxelCubeGeneration.h>
#include <spdlog/spdlog.h>
#include <SOAPSettings.h>
#include <AsciiArt.h>
#include <fstream>
#include <IPosition.h>

using namespace YAML;

void validateClusteringSettings(const YAML::Node &inputYaml) {

    auto clusteringNode = inputYaml["Clustering"];

    for (auto node : clusteringNode) {
        auto methodName = node.first.as<std::string>();

        switch (IClusterer::typeFromString(methodName)) {
            case IClusterer::Type::IdentityClusterer: {
                IdentityClusterer::settings
                        = Settings::IdentityClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::DistanceClusterer: {
                DistanceClusterer::settings
                        = Settings::DistanceClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::DensityBasedClusterer: {
                DensityBasedClusterer::settings
                        = Settings::DensityBasedClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::ReferencePositionsClusterer: {
                ReferencePositionsClusterer::settings
                        = Settings::ReferencePositionsClusterer(clusteringNode);
                break;
            }
            case IClusterer::Type::SOAPClusterer: {
                SOAPClusterer::settings
                        = Settings::SOAPClusterer(clusteringNode);

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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.yml" << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];
    std::string::size_type idx;

    idx = inputFilename.rfind('.');

    std::string filenameWithoutExtension;
    if(idx != std::string::npos)
        filenameWithoutExtension = inputFilename.substr(idx);

    YAML::Node inputYaml = YAML::LoadFile(argv[1]);
    YAML::Emitter emitter;
    emitter << inputYaml;

    YAML::Emitter outputYaml;
    spdlog::info("Executable: {}", argv[0]);
    spdlog::info("Validating input from file {}:\n{}...", inputFilename, emitter.c_str());
    validateInput(inputYaml);

    // Apply settings from inputYaml
    MaximaProcessing::settings = Settings::MaximaProcessing(inputYaml);
    unsigned samplesToAnalyze = MaximaProcessing::settings.samplesToAnalyze();
    if (samplesToAnalyze == 0)
        spdlog::info("Analyzing all samples.");
    else if(samplesToAnalyze < std::numeric_limits<unsigned>::max())
        spdlog::info("Analyzing {} samples.", samplesToAnalyze);

    // if binary file basename is not specified, use the name of  the input file
    if(!inputYaml["MaximaProcessing"][MaximaProcessing::settings.binaryFileBasename.name()]) {
        idx = inputFilename.rfind('.');
        if (idx != std::string::npos)
            MaximaProcessing::settings.binaryFileBasename = inputFilename.substr(0,idx);
    }
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
            case IClusterer::Type::IdentityClusterer: {
                auto &settings = IdentityClusterer::settings;

                settings = Settings::IdentityClusterer(node.second);

                IdentityClusterer identityClusterer(samples);
                if (!node.second[settings.identityValueIncrement.name()])
                    settings.identityValueIncrement = valueStandardError * 1e-4;

                identityClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::DistanceClusterer: {
                auto &settings = DistanceClusterer::settings;

                settings = Settings::DistanceClusterer(node.second);

                DistanceClusterer distanceClusterer(samples);
                if (!node.second[settings.similarityValueIncrement.name()])
                    settings.similarityValueIncrement = valueStandardError * 1e-2;
                distanceClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::DensityBasedClusterer: {
                auto &settings = DensityBasedClusterer::settings;

                settings = Settings::DensityBasedClusterer(node.second);

                DensityBasedClusterer densityBasedClusterer(samples);
                densityBasedClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::ReferencePositionsClusterer: {
                auto &settings = ReferencePositionsClusterer::settings;

                settings = Settings::ReferencePositionsClusterer(node.second);

                std::vector<Eigen::Vector3d> positions;

                auto positionNodes = node.second["positions"];
                for (const auto &positionNode : positionNodes){
                    positions.emplace_back(YAML::decodePosition(positionNode, atoms));
                }

                spdlog::info("Using the following positions:");
                for (const auto &position : positions){
                    spdlog::info("{} {} {}", position[0], position[1], position[2]);
                }

                ReferencePositionsClusterer ReferencePositionsClusterer(samples, atoms, positions);
                ReferencePositionsClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IClusterer::Type::SOAPClusterer: {
                auto &settings = SOAPClusterer::settings;

                settings = Settings::SOAPClusterer(node.second);

                auto soapSettings = node.second["SOAP"];
                if (soapSettings[SOAP::General::settings.name()])
                    SOAP::General::settings = Settings::SOAP::General(soapSettings);
                if (soapSettings[SOAP::Radial::settings.name()])
                    SOAP::Radial::settings = Settings::SOAP::Radial(soapSettings);
                if (soapSettings[SOAP::Angular::settings.name()])
                    SOAP::Angular::settings = Settings::SOAP::Angular(soapSettings);
                if (soapSettings[SOAP::Cutoff::settings.name()])
                    SOAP::Cutoff::settings = Settings::SOAP::Cutoff(soapSettings);

                SOAPClusterer sOAPClusterer(atoms, samples);
                sOAPClusterer.cluster(maxima);

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
    spdlog::info("Writing results into file \"{}\"", resultsFilename);

    std::ofstream yamlFile(resultsFilename);
    yamlFile << outputYaml.c_str();
    yamlFile.close();

    spdlog::info("Done! Bye bye.");

    return 0;

    /* TODO
     * - refactor best match clusterer
     * - make single value statistics class
     * - test naive standard deviation
     * - test choice of function value increment
     * - split identity sort into batches that can be compared in parallel using OpenMP
     * */
};
