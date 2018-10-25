//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <GlobalClusterSorter.h>
#include <EnergyCalculator.h>

#include <MoleculeWidget.h>
#include <ElectronicWaveFunction.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <QApplication>
#include <algorithm>
#include <utility>

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    //auto atoms = ElectronicWaveFunction::getInstance(std::string("Acetone-em.wf")).getAtomsVector();
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("Acetone.bin",100);
    auto atoms = reader.getAtoms();

    console->info("number of inital refs {}",globallyIdenticalMaxima.size());

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

    //Statistics
    energyCalculator.calculateStatistics(globallyClusteredMaxima);
    std::ofstream yamlFile("energies.yml");
    yamlFile << energyCalculator.getYamlDocumentString();
    yamlFile.close();

    // Visuals
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();
    AtomsVector3D(root, atoms);

    for (auto i : globallyClusteredMaxima[0]){
        ElectronsVector3D(root, atoms, i.representativeReference().maximum(), false);
    }
    return QApplication::exec();

};
