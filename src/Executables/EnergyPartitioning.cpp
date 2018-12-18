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

void checkBasename(char *arg, std::string &basename) {
    basename = arg;
    try {
        if (basename.find('.') != std::string::npos)
            throw std::invalid_argument("Basename '" + basename + "'should not have a file ending.");
    }
    catch (std::invalid_argument &e){
        std::cout << e.what() << std::endl;
        abort();
    }
};

void checkNumberOfSamples(char *arg, int &numberOfSamples) {
    numberOfSamples = std::atoi(arg);
    try {
        if (numberOfSamples < 1)
            throw std::invalid_argument("NumberOfSamples is '" + std::to_string(numberOfSamples) + "'but should be positive.");
    }
    catch (std::invalid_argument &e){
        std::cout << e.what() << std::endl;
        abort();
    }
};

bool handleCommandlineArguments(int argc, char *const *argv,
                                std::string &basename,
                                int &numberOfSamples) {
    if (argc < 2) {
        std::cout << "Usage: \n"
                  << "Argument 1: basename\n"
                  << "Argument 2: number of samples (optional)"
                  << std::endl;
        std::cout << "raw 100" << std::endl;
        return false;
    } else if (argc == 2) {
        checkBasename(argv[1], basename);
        return true;
    } else if (argc == 3) {
        checkBasename(argv[1], basename);
        checkNumberOfSamples(argv[2], numberOfSamples);
        return true;
    } else {
        throw std::invalid_argument("Too many arguments");
    }
}


int main(int argc, char *argv[]) {
    std::string basename;
    int numberOfSamples = std::numeric_limits<int>::max();

    if (basename.empty()) {
        bool inputArgumentsFoundQ =
                handleCommandlineArguments(argc, argv, basename, numberOfSamples);
        if (!inputArgumentsFoundQ) return 1;
    }

    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read(basename, size_t(numberOfSamples));
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


    double identicalDistThresh = 0.01;
    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.sort();
    console->info("number of elements after identity sort {}",globallyIdenticalMaxima.size());

    double similarDistThresh = 0.2;
    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(samples,globallyIdenticalMaxima, globallySimilarMaxima,similarDistThresh);
    globalSimilaritySorter.sort();
    console->info("number of elements after similarity sort {}",globallySimilarMaxima.size());

    std::vector<std::vector<SimilarReferences>> globallyClusteredMaxima;
    GlobalClusterSorter globalClusterSorter(samples, globallySimilarMaxima, globallyClusteredMaxima, similarDistThresh);
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
     * - load settings from .yaml file
     * */
};
