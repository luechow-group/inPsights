//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <GlobalClusterSorter.h>
#include <EnergyCalculator.h>
#include <GeneralStatistics.h>
#include <algorithm>
#include <utility>
#include <EnergyPartitioningSettings.h>
#include <VoxelCubeGeneration.h>
#include <SurfaceDataGenerator.h>

using namespace YAML;
using namespace Logger;

int main(int argc, char *argv[]) {

    std::string inputFilename = argv[1];
    YAML::Node inputYaml = YAML::LoadFile(argv[1]);
    YAML::Emitter emitter;
    emitter << inputYaml;

    YAML::Emitter outputYaml;
    console->info("Executable: {}", argv[0]);
    console->info("Input file {}:\n{}", inputFilename, emitter.c_str());

    // Apply settings from inputYaml
    Settings::EnergyPartitioning settings(inputYaml);
    GlobalIdentitySorter::settings = Settings::GlobalIdentitySorter(inputYaml);
    GlobalSimilaritySorter::settings = Settings::GlobalSimilaritySorter(inputYaml);
    GlobalClusterSorter::settings = Settings::GlobalClusterSorter(inputYaml);
    VoxelCubeGeneration::settings = Settings::VoxelCubeGeneration(inputYaml);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima, samples);
    reader.read(settings.binaryFileBasename.get(), settings.samplesToAnalyze.get());
    auto atoms = reader.getAtoms();

    console->info("number of inital refs {}", globallyIdenticalMaxima.size());
    auto results = GeneralStatistics::calculate(globallyIdenticalMaxima, samples, atoms);

    EnergyCalculator energyCalculator(outputYaml, samples,atoms);

    auto valueStandardError = results.valueStats_.standardError()(0,0);

    if(settings.identitySearch.get()) {
        console->info("Start identity search");
        GlobalIdentitySorter globalIdentiySorter(globallyIdenticalMaxima,samples);
        if(!inputYaml["GlobalIdentitySorter"]["identityValueIncrement"])
            GlobalIdentitySorter::settings.identityValueIncrement = valueStandardError*1e-4;
        globalIdentiySorter.sort();
        console->info("number of elements after identity sort {}", globallyIdenticalMaxima.size());
    }

    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(samples,globallyIdenticalMaxima, globallySimilarMaxima);
    if(!inputYaml["GlobalSimilaritySorter"]["similarityValueIncrement"])
        GlobalSimilaritySorter::settings.similarityValueIncrement = valueStandardError*1e-2;
    globalSimilaritySorter.sort();
    console->info("number of elements after similarity sort {}", globallySimilarMaxima.size());

    std::vector<std::vector<SimilarReferences>> globallyClusteredMaxima;
    GlobalClusterSorter globalClusterSorter(
            samples,
            globallySimilarMaxima,
            globallyClusteredMaxima);
    globalClusterSorter.sort();
    console->info("number of elements after cluster sort {}", globallyClusteredMaxima.size());

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
    console->info("Calculating statistics...");
    energyCalculator.calculateStatistics(globallyClusteredMaxima);
    //outputYaml << Key << "SurfaceData" << Value << surfaceData;

    if(VoxelCubeGeneration::settings.generateVoxelCubesQ.get()) {
        auto voxelCubes = VoxelCubeGeneration::fromCluster(globallyClusteredMaxima[0], samples);
        outputYaml << Key << "VoxelData" << Value << voxelCubes;
    }

    outputYaml << EndMap << EndDoc;

    std::string resultsFilename = settings.binaryFileBasename.get() + ".yml";
    console->info("Writing results into file \"{}\"", resultsFilename);


    std::ofstream yamlFile(resultsFilename);
    yamlFile << outputYaml.c_str();
    yamlFile.close();

    console->info("Done! Bye bye.");

    return 0;

    /*TODO
     * - make similar ref a nested structure
     *  - has a first/representative structure
     *  - has an averagedStructure
     *  - similarity attribute (enum) - spatially|permutationally|rotationally + similar|identical,
     * - use averaged structure for permutation sort, store best match permutation to add energies
     * - use spin correlations to select motifs
     *  - calculate inter- and intra- motif energies
     *   - what happens with the cores?
     * - export voxel for all electrons for all clusters (external file?)
     *  - center voxel cube on averaged electron pos
     * - refactor names
     *  - rename EnergyPartitioning to MaximaProcessing
     * - settings classes
     *   - fix memory leak after destruction of properties observed with the intel compiler
     *   - voxel settings
     *   - expansion settings structural similarity
     * - make single value statistics class
     * - F2 cluster includes other references that shouldn't be there
     * - F2 spin correlations seem to be not symmetrical
     * - test naive std
     * - choice of function value increment
     * - validate that ring-like clusters are ordered correctly
     * - use global similarity for permutation sorting
     * - split identity sort into batches that can be compared in parallel using OpenMP
     * - improve spinSpecificHungarian (low priority)
     * */
};
