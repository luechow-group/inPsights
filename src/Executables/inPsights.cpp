//
// Created by heuer on 15.11.18.
//

#include <iostream>
#include <MoleculeWidget.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <QApplication>
#include <yaml-cpp/yaml.h>
#include <ParticlesVector.h>

bool handleCommandlineArguments(int argc, char *const *argv,
                                std::string &resultFilename,
                                size_t &clusterId) {
    if (argc < 3) {
        std::cout << "Usage: \n"
                  << "Argument 1: result file\n"
                  << "Argument 2: cluster id"
                  << std::endl;
        std::cout << "raw.yml 100" << std::endl;
        return false;
    } else if (argc == 3) {
        clusterId = std::atol(argv[2]);
        return true;
    } else {
        throw std::invalid_argument("Too many arguments");
    }
}

int main(int argc, char *argv[]) {
    std::string resultFilename;
    size_t clusterId = 0;

    if (resultFilename.empty()) {
        bool inputArgumentsFoundQ = handleCommandlineArguments(argc, argv, resultFilename, clusterId);
        if (!inputArgumentsFoundQ) return 1;
    }

    YAML::Node doc = YAML::LoadFile("raw.yml");
    auto atoms = doc["Atoms"].as<AtomsVector>();
    auto Vnn = doc["Vnn"];

    auto nClusters =doc["Clusters"].size();
    assert(clusterId < nClusters && "The cluster id must be smaller than the total number of clusters.");
    std::cout << "analyzing cluster with id " << clusterId << " out of [0," << nClusters-1 << "]" << std::endl;

    auto cluster = doc["Clusters"][clusterId];
    auto electronsVectorCollection = cluster["Structures"].as<std::vector<ElectronsVector>>();
    auto spinCorrelations = cluster["SpinCorrelations"];
    auto Te = cluster["Te"];
    auto Vee = cluster["Vee"];
    auto Ven = cluster["Ven"];

    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();
    AtomsVector3D(root, atoms);
    ElectronsVector3D(root, electronsVectorCollection[0]);
    //for (const auto & i : electronsVectorCollection)
    //    ElectronsVector3D(root, i);

    auto nElectrons = Te.size();
    std::vector<double> Ve(nElectrons);
    std::vector<double> VeErr(nElectrons);

    for (int i = 0; i < nElectrons; ++i) {
        Ve[i] = Te[i][0].as<double>();
        VeErr[i] = Te[i][1].as<double>();

        for (int k = 0; k < atoms.numberOfEntities(); ++k)
            Ve[i] += Ven[i][k][0].as<double>();

        for (int j = i+1; j < nElectrons; ++j)
            Ve[i] += Vee[i][j][0].as<double>();

    }
    for (auto e : Ve)
        std::cout << e << std::endl;

    return QApplication::exec();

    /* TODO
     * YAML
     * proper conversion of statistic nodes (templatized)
     *
     * METHOD
     * - spatial permutations (by value range or struct sim)
     *
     * GUI
     * - correlation in plots via button
     *  - dotted and dashed lines
     * - same spin connections via button
     * - select clusters
     * - show energies
     * - show eigenvectors (.wf needed)
     */
};