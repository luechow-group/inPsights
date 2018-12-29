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
#include <experimental/filesystem>
#include <GlobalSortSettings.h>

using namespace YAML;
using namespace Logger;

int main(int argc, char *argv[]) {
    YAML::Node doc = YAML::LoadFile(argv[1]);
    YAML::Emitter emitter;
    emitter << doc;
    console->info("Executable: {}", argv[0]);
    console->info("Input file {}:\n{}", argv[1], emitter.c_str());

    Settings::GlobalSort settings(doc);
    auto basename = doc["binaryFileBasename"].as<std::string>();

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima, samples);
    reader.read(basename, settings.samplesToAnalyze.get());
    auto atoms = reader.getAtoms();

    console->info("number of inital refs {}", globallyIdenticalMaxima.size());
    auto results = GeneralStatistics::calculate(globallyIdenticalMaxima, samples, atoms);

    YAML::Emitter out;
    out << BeginDoc << BeginMap
        << Key << "Atoms" << Value << atoms << Comment("[a0]")
        << Key << "NSamples" << Value << samples.size()
        << Key << "OverallResults" << results;

    EnergyCalculator energyCalculator(out, samples,atoms);

    auto valueStandardError = results.valueStats_.standardError()(0,0);

    if(settings.identitySearch.get()) {
        console->info("Start identity search");
        GlobalIdentiySorter globalIdentiySorter(
                globallyIdenticalMaxima,
                samples,
                settings.identityRadius.get(),
                valueStandardError*1e-4);
        globalIdentiySorter.sort();
        console->info("number of elements after identity sort {}",globallyIdenticalMaxima.size());
    }

    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(samples,
            globallyIdenticalMaxima,
            globallySimilarMaxima,
            settings.similarityRadius.get(), valueStandardError*1e-2);
    globalSimilaritySorter.sort();
    console->info("number of elements after similarity sort {}", globallySimilarMaxima.size());

    std::vector<std::vector<SimilarReferences>> globallyClusteredMaxima;
    GlobalClusterSorter globalClusterSorter(
            samples,
            globallySimilarMaxima,
            globallyClusteredMaxima,
            settings.similarityRadius.get());
    globalClusterSorter.sort();
    console->info("number of elements after cluster sort {}", globallyClusteredMaxima.size());

    // Permutation sort
    /*std::vector<std::vector<SimilarReferences>> globallyPermutationallyInvariantClusteredMaxima;
    GlobalPermutationSorter globalPermutationSorter(atoms, samples, globallyClusteredMaxima, globallyPermutationallyInvariantClusteredMaxima);
    globalPermutationSorter.sort();
    */

    console->info("Calculating statistics...");
    energyCalculator.calculateStatistics(globallyClusteredMaxima);

    std::string resultsFilename = basename + ".yml";
    console->info("Writing results into file \"{}\"",resultsFilename);
    out << EndDoc << EndMap;
    std::ofstream yamlFile(resultsFilename);
    yamlFile << out.c_str();
    yamlFile.close();
    console->info("Done! Bye bye.");

    return 0;

    /*TODO
     * - make single value statistics class
     * - refactor names
     * - F2 cluster includes other references that shouldn't be there
     * - test naive std
     * - choice of function value increment
     * - validate that ring-like clusters are ordered correctly
     * - use global similarity for permutation sorting
     * - split identity sort into batches that can be compared in parallel using OpenMP
     * - improve spinSpecificHungarian
     * */
};
