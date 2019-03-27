//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <GlobalClusterSorter.h>
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
    Settings::MaximaProcessing settings(inputYaml);
    GlobalIdentitySorter::settings = Settings::GlobalIdentitySorter(inputYaml);
    GlobalSimilaritySorter::settings = Settings::GlobalSimilaritySorter(inputYaml);
    GlobalClusterSorter::settings = Settings::GlobalClusterSorter(inputYaml);
    VoxelCubeGeneration::settings = Settings::VoxelCubeGeneration(inputYaml);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima, samples);
    reader.read(settings.binaryFileBasename(), settings.samplesToAnalyze());
    auto atoms = reader.getAtoms();

    spdlog::info("number of inital refs {}", globallyIdenticalMaxima.size());
    auto results = GeneralStatistics::calculate(globallyIdenticalMaxima, samples, atoms);

    auto valueStandardError = results.valueStats_.standardError()(0,0);

    if(settings.identitySearch()) {
        spdlog::info("Start identity search");
        GlobalIdentitySorter globalIdentiySorter(globallyIdenticalMaxima,samples);
        if(!inputYaml["GlobalIdentitySorter"]["identityValueIncrement"])
            GlobalIdentitySorter::settings.identityValueIncrement = valueStandardError*1e-4;
        globalIdentiySorter.sort();
        spdlog::info("number of elements after identity sort {}", globallyIdenticalMaxima.size());
    }

    std::vector<SimilarReferences> globallySimilarMaxima;

    GlobalSimilaritySorter globalSimilaritySorter(samples,globallyIdenticalMaxima, globallySimilarMaxima);
    if(!inputYaml["GlobalSimilaritySorter"]["similarityValueIncrement"])
        GlobalSimilaritySorter::settings.similarityValueIncrement = valueStandardError*1e-2;
    globalSimilaritySorter.sort();
    spdlog::info("number of elements after similarity sort {}", globallySimilarMaxima.size());

    std::vector<std::vector<SimilarReferences>> globallyClusteredMaxima;
    GlobalClusterSorter globalClusterSorter(
            samples,
            globallySimilarMaxima,
            globallyClusteredMaxima);
    globalClusterSorter.sort();
    spdlog::info("number of elements after cluster sort {}", globallyClusteredMaxima.size());

    // Permutation sort
    /*std::vector<std::vector<SimilarReferences>> globallyPermutationallyInvariantClusteredMaxima;
    GlobalPermutationSorter globalPermutationSorter(atoms, samples, globallyClusteredMaxima, globallyPermutationallyInvariantClusteredMaxima);
    globalPermutationSorter.sort();
    */

    // write used settings
    YAML::Node usedSettings;
    settings.appendToNode(usedSettings);
    GlobalIdentitySorter::settings.appendToNode(usedSettings);
    GlobalSimilaritySorter::settings.appendToNode(usedSettings);
    GlobalClusterSorter::settings.appendToNode(usedSettings);
    VoxelCubeGeneration::settings.appendToNode(usedSettings);
    outputYaml << BeginDoc << Comment("used settings") << usedSettings << EndDoc;

    // write results
    outputYaml << BeginDoc << BeginMap
               << Key << "Atoms" << Value << atoms << Comment("[a0]")
               << Key << "NSamples" << Value << samples.size()
               << Key << "OverallResults" << Value << results;
    spdlog::info("Calculating statistics...");

    MaximaProcessor energyCalculator(outputYaml, samples,atoms);
    energyCalculator.calculateStatistics(globallyClusteredMaxima);

    outputYaml << EndMap << EndDoc;

    std::string resultsFilename = settings.binaryFileBasename() + ".yml";
    spdlog::info("Writing results into file \"{}\"", resultsFilename);

    std::ofstream yamlFile(resultsFilename);
    yamlFile << outputYaml.c_str();
    yamlFile.close();

    spdlog::info("Done! Bye bye.");

    return 0;

    /*TODO
     * - make similar ref a nested structure
     *  - has a first/representative structure
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
