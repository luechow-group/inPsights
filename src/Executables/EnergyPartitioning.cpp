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

namespace fs = std::experimental::filesystem;

bool handleCommandlineArguments(int argc, char *const *argv, std::string &fileName) {
    if (argc < 2) {
        std::cout << "Usage: \n"
                  << "Argument 1: filename\n"
                  << std::endl;
        std::cout << "epart.yml" << std::endl;
        return false;
    } else if (argc == 2) {
        fileName = argv[1];
        auto extension = fs::path(fileName).extension().string();
        std::vector<std::string> extensions = {".yml", ".yaml", ".json"};
        return std::any_of(extensions.begin(), extensions.end(), [extension](std::string i){return i == extension;});
    } else {
        throw std::invalid_argument("Too many arguments");
    }
}

using namespace YAML;

int main(int argc, char *argv[]) {
    std::string fileName;
    if (fileName.empty()) {
        bool inputArgumentsFoundQ = handleCommandlineArguments(argc, argv, fileName);
        if (!inputArgumentsFoundQ)
            return 1;
    }

    Logger::initialize();
    auto console = Logger::get();

    YAML::Node doc = YAML::LoadFile(fileName);

    //print input file
    YAML::Emitter emitter;
    emitter << doc;
    console->info("Executable: {}", argv[0]);
    console->info("Input file {}:\n{}", fileName, emitter.c_str());

    size_t numberOfSamples;
    if(doc["settings"]["samplesToAnalyze"])
        try {
            numberOfSamples = doc["settings"]["samplesToAnalyze"].as<size_t>();
            console->info("Analyzing {} samples.", numberOfSamples);
        } catch ( const YAML::BadConversion&  e) {
            auto string = doc["settings"]["samplesToAnalyze"].as<std::string>();
            if(string == "all" || string == "All" || string == "ALL"){
                numberOfSamples = std::numeric_limits<size_t>::max();
                console->info("Analyzing all {} samples.", numberOfSamples);
            } else {
                console->error("{0} \n"
                               "The key samplesToAnalyze is neither a positive integer nor \"all\". Abort.", e.what());
                return 1;
            }
        }
    else {
        numberOfSamples = std::numeric_limits<size_t>::max();
        console->info("Analyzing all {} samples.", numberOfSamples);
    }
    auto basename = doc["binaryFileBasename"].as<std::string>();

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima, samples);
    reader.read(basename, numberOfSamples);
    auto atoms = reader.getAtoms();

    console->info("number of inital refs {}", globallyIdenticalMaxima.size());
    auto results = GeneralStatistics::calculate(globallyIdenticalMaxima, samples, atoms);

    YAML::Emitter out;
    out << BeginDoc << BeginMap
        << Key << "Atoms" << Value << atoms << Comment("[a0]")
        << Key << "NSamples" << Value << samples.size()
        << Key << "OverallResults" << results;

    EnergyCalculator energyCalculator(out, samples,atoms);

    auto identityRadius = doc["settings"]["identityRadius"].as<double>();
    auto similarityRadius = doc["settings"]["similarityRadius"].as<double>();
    auto valueStandardError = results.valueStats_.standardError()(0,0);

    auto identitySearchOption = doc["settings"]["skipIdentitySearch"];
    if(identitySearchOption && identitySearchOption.as<bool>()) {
        console->info("Start identity search");
        GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identityRadius, valueStandardError*1e-4);
        globalIdentiySorter.sort();
        console->info("number of elements after identity sort {}",globallyIdenticalMaxima.size());
    }

    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(samples,
            globallyIdenticalMaxima,
            globallySimilarMaxima,
            similarityRadius, valueStandardError*1e-2);
    globalSimilaritySorter.sort();
    console->info("number of elements after similarity sort {}", globallySimilarMaxima.size());

    std::vector<std::vector<SimilarReferences>> globallyClusteredMaxima;
    GlobalClusterSorter globalClusterSorter(samples, globallySimilarMaxima, globallyClusteredMaxima, similarityRadius);
    globalClusterSorter.sort();
    console->info("number of elements after cluster sort {}", globallyClusteredMaxima.size());

    // Permutation sort
    /*std::vector<std::vector<SimilarReferences>> globallyPermutationallyInvariantClusteredMaxima;
    GlobalPermutationSorter globalPermutationSorter(atoms, samples, globallyClusteredMaxima, globallyPermutationallyInvariantClusteredMaxima);
    globalPermutationSorter.sort();
    */

    //Statistics
    energyCalculator.calculateStatistics(globallyClusteredMaxima);

    out << EndDoc << EndMap;
    std::ofstream yamlFile(basename + ".yml");
    yamlFile << out.c_str();
    yamlFile.close();
    
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
