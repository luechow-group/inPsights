//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <GlobalClusterSorter.h>
#include <GlobalPermutationSorter.h>
#include <EnergyCalculator.h>
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


int main(int argc, char *argv[]) {
    std::string fileName;
    if (fileName.empty()) {
        bool inputArgumentsFoundQ = handleCommandlineArguments(argc, argv, fileName);
        if (!inputArgumentsFoundQ)
            return 1;
    }

    YAML::Node doc = YAML::LoadFile(fileName);
    auto numberOfSamples = doc["settings"]["samplesToAnalyze"].as<size_t >();
    auto basename = doc["basename"].as<std::string>();


    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read(basename, numberOfSamples);
    auto atoms = reader.getAtoms();

    console->info("number of inital refs {}", globallyIdenticalMaxima.size());

    EnergyCalculator energyCalculator(samples,atoms);
    auto totalEnergies = energyCalculator.calculateTotalEnergies();

    console->info("Te= {}, Vee = {}, Ven = {}, Vnn = {}, Eges = {}",
            totalEnergies.Te,
            totalEnergies.Vee,
            totalEnergies.Ven,
            totalEnergies.Vnn,
            totalEnergies.totalEnergy());

    auto identityRadius = doc["settings"]["identityRadius"].as<double>();
    auto similarityRadius = doc["settings"]["similarityRadius"].as<double>();

    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identityRadius);
    globalIdentiySorter.sort();
    console->info("number of elements after identity sort {}",globallyIdenticalMaxima.size());

    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(samples,globallyIdenticalMaxima, globallySimilarMaxima,similarityRadius);
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
    std::ofstream yamlFile(basename + ".yml");
    yamlFile << energyCalculator.getYamlDocumentString();
    yamlFile.close();
    
    return 0;

    /*TODO
     * - F2 cluster includes other references that shouldn't be there
     * - test naive std
     * - choice of function value increment
     * - validate that ring-like clusters are ordered correctly
     * - use global similarity for permutation sorting
     * */
};
