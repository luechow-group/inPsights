// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <RawDataReader.h>
#include <Cluster.h>
#include <ElementInfo.h>
#include <IdentityClusterer.h>
#include <PreClusterer.h>
#include <DensityBasedClusterer.h>
#include <SOAPClusterer.h>
#include <ReferencePositionsClusterer.h>
#include <ElectronSelection.h>
#include <MaximaProcessor.h>
#include <GeneralStatistics.h>
#include <algorithm>
#include <MaximaProcessingSettings.h>
#include "CameraSettings.h"
#include <VoxelCubeGeneration.h>
#include <VoxelCubeOverlapCalculation.h>
#include <spdlog/spdlog.h>
#include <SOAPSettings.h>
#include <Logo.h>
#include <fstream>
#include <Metrics.h>
#include <MolecularSelection.h>

using namespace YAML;

void validateClusteringSettings(const YAML::Node &inputYaml) {

    auto clusteringNode = inputYaml["Clustering"];

    for (const auto& node : clusteringNode) {
        auto methodName = node.first.as<std::string>();

        switch (IProcess::typeFromString(methodName)) {
            case IProcess::ProcessType::IdentityClusterer: {
                IdentityClusterer::settings = Settings::IdentityClusterer(clusteringNode);
                break;
            }
            case IProcess::ProcessType::DistanceClusterer: {
                PreClusterer::settings = Settings::PreClusterer(clusteringNode);
                break;
            }
            case IProcess::ProcessType::DensityBasedClusterer: {
                DensityBasedClusterer::settings = Settings::DensityBasedClusterer(clusteringNode);
                break;
            }
            case IProcess::ProcessType::ReferencePositionsClusterer: {
                ReferencePositionsClusterer::settings = Settings::ReferencePositionsClusterer(clusteringNode);
                break;
            }
            case IProcess::ProcessType::SOAPClusterer: {
                SOAPClusterer::settings = Settings::SOAPClusterer(clusteringNode);

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
    spdlog::info("Version: {}", inPsights::version());
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
            MaximaProcessing::settings.binaryFileBasename = inputFilename.substr(0, idx);
    }
        // Read maxima and samples
    Cluster maxima;
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

    maxima.sortAll();
    for (const auto& node : clusteringNode) {
        auto methodName = node.first.as<std::string>();
        spdlog::info("Starting {}...", methodName);

        if (usedClusteringSettings[methodName])
            spdlog::warn("Method \"{}\" is being applied multiple times! Overwriting old settings...", methodName);

        switch (IProcess::typeFromString(methodName)) {
            case IProcess::ProcessType::IdentityClusterer: {
                auto &settings = IdentityClusterer::settings;

                settings = Settings::IdentityClusterer(node.second);

                IdentityClusterer identityClusterer(samples);
                if (!node.second[settings.valueIncrement.name()])
                    settings.valueIncrement = valueStandardError * 1e-4;

                identityClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IProcess::ProcessType::DistanceClusterer: {
                auto &settings = PreClusterer::settings;

                settings = Settings::PreClusterer(node.second);

                PreClusterer distanceClusterer(samples);
                if (!node.second[settings.valueIncrement.name()])
                    settings.valueIncrement = valueStandardError * 1e-2;
                distanceClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IProcess::ProcessType::DensityBasedClusterer: {
                auto &settings = DensityBasedClusterer::settings;

                settings = Settings::DensityBasedClusterer(node.second);

                if(settings.local()) {
                    auto electronSelectionSettings = node.second[VARNAME(ElectronSelection)];
                    if (electronSelectionSettings)
                        ElectronSelection::settings = Settings::ElectronSelection(electronSelectionSettings, atoms);
                }

                DensityBasedClusterer densityBasedClusterer(samples);
                densityBasedClusterer.cluster(maxima);

                if(settings.local()) {
                    ElectronSelection::settings.appendToNode(usedClusteringSettings); // TODO FIX USED SETTINGS APPEND
                }

                settings.appendToNode(usedClusteringSettings);
                break;
            }
            case IProcess::ProcessType::ReferencePositionsClusterer: {
                auto &settings = ReferencePositionsClusterer::settings;

                settings = Settings::ReferencePositionsClusterer(node.second);

                if(settings.local()) {
                    auto electronSelectionSettings = node.second[VARNAME(ElectronSelection)];
                    if (electronSelectionSettings)
                        ElectronSelection::settings = Settings::ElectronSelection(electronSelectionSettings, atoms);
                }
                ReferencePositionsClusterer ReferencePositionsClusterer(samples);
                ReferencePositionsClusterer.cluster(maxima);

                settings.appendToNode(usedClusteringSettings);

                if(settings.local()) {
                    ElectronSelection::settings.appendToNode(usedClusteringSettings); // TODO FIX USED SETTINGS APPEND
                }

                break;
            }
            case IProcess::ProcessType::SOAPClusterer: {
                auto &settings = SOAPClusterer::settings;

                settings = Settings::SOAPClusterer(node.second);

                auto soapSettings = node.second["SOAP"]; //TODO what if node "SOAP" is not present?
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

    VoxelCubeOverlapCalculation::settings = Settings::VoxelCubeOverlapCalculation(inputYaml);
    VoxelCubeOverlapCalculation::settings.appendToNode(usedSettings);

    outputYaml << BeginDoc
    << Comment("input from \"" + inputFilename + "\"") << inputYaml << EndDoc;
    outputYaml << BeginDoc
               << BeginMap
               << Comment(inPsights::logo);

    YAML::Node cameraNode;
    Camera::settings.appendToNode(cameraNode);

    outputYaml << Key << "Camera" << Value << cameraNode["Camera"] << Comment("[%,°,°,°]");

    outputYaml << Key << "Atoms" << Value << atoms << Comment("[a0]")
               << Key << "Rnn" << Value << Metrics::positionalDistances(atoms.positionsVector()) << Comment("[a0]")
               << Key << "NSamples" << Value << samples.size()
               << Key << "OverallResults" << Value << results;

    spdlog::info("Calculating statistics...");

    MaximaProcessor maximaProcessor(outputYaml, samples, atoms);

    // TODO REFACTOR
    std::vector<std::vector<std::vector<size_t>>> nucleiMergeLists;
    if(inputYaml["MotifMerge"]) {
        auto motifMergeNode = inputYaml["MotifMerge"];
        for (YAML::iterator it = motifMergeNode.begin(); it != motifMergeNode.end(); ++it) {
            const YAML::Node &nucleiMergeListNode = *it;

            auto nucleiMergeList = it->as<std::vector<std::vector<size_t>>>();
            nucleiMergeLists.emplace_back(nucleiMergeList);
        }
    }

    std::vector<DynamicMolecularSelection> selections;
    if(inputYaml["EnergySelections"]) {
        auto energySelectionNode = inputYaml["EnergySelections"];
        for (YAML::iterator it = energySelectionNode.begin(); it != energySelectionNode.end(); ++it) {
            const YAML::Node &node = *it;

            auto nucleiIndices = node["Nuclei"].as<ParticleIndices>();
            auto nearestElectronsSettings = Settings::ElectronSelection(node["ElectronSelection"], atoms);

            selections.emplace_back(DynamicMolecularSelection(nearestElectronsSettings, nucleiIndices));
        }
    }


    // local particle energies
    std::vector<size_t> nucleiIndices(0);
    if(inputYaml["LocalParticleEnergies"]) {
        auto collectiveEnergyNode = inputYaml["LocalParticleEnergies"];
        nucleiIndices = collectiveEnergyNode["nuclei"].as<std::vector<size_t>>();
    }

    maximaProcessor.calculateStatistics(maxima, nucleiMergeLists, nucleiIndices, selections);

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
};
