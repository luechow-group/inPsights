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
#include <Statistics.h>
#include <InPsightsWidget.h>


void singleElectronEnergies(const YAML::Node &cluster) {
    auto Te = cluster["Te"];
    auto Vee = cluster["Vee"];
    auto Ven = cluster["Ven"];
    auto nElectrons = Te.size();
    auto nAtoms = Ven.begin()->second.size();
    std::vector<double> Ve(nElectrons);
    std::vector<double> VeErr(nElectrons);

    for (int i = 0; i < nElectrons; ++i) {
        Ve[i] = Te[i][0].as<double>();
        VeErr[i] = Te[i][1].as<double>();

        for (int k = 0; k < nAtoms; ++k){
            Ve[i] += Ven[i][k][0].as<double>();
            VeErr[i] += Ven[i][k][1].as<double>();
        }

        for (int j = i+1; j < nElectrons; ++j){
            Ve[i] += 0.5*Vee[i][j][0].as<double>();
            VeErr[i] += 0.5*Vee[i][j][1].as<double>();
        }
        //TODO FEHLERADDITION
        std::cout << "Ee = " << Ve[i] << "+/-" << VeErr[i] << "  (ERROR is wrong -> add error propagation)" << std::endl;
    }
};

bool handleCommandlineArguments(int argc, char *const *argv, std::string &resultFilename) {
    if (argc < 2) {
        std::cout << "Usage: \n"
                  << "Argument 1: result file\n"
                  << std::endl;
        std::cout << "raw.yml" << std::endl;
        return false;
    } else if (argc == 2) {
        resultFilename = argv[1];
        return true;
    } else {
        throw std::invalid_argument("Too many arguments");
    }
}


int main(int argc, char *argv[]) {
    std::string resultFilename;
    if (resultFilename.empty()) {
        bool inputArgumentsFoundQ = handleCommandlineArguments(argc, argv, resultFilename);
        if (!inputArgumentsFoundQ) return 1;
    }

    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    YAML::Node doc = YAML::LoadFile("raw.yml");
    size_t clusterId = 0;
    auto cluster = doc["Clusters"][clusterId];

    auto atoms = doc["Atoms"].as<AtomsVector>();
    auto electronsVectorCollection = cluster["Structures"].as<std::vector<ElectronsVector>>();
    auto spinCorrelations = cluster["SpinCorrelations"];

    auto inPsightsWidget = new InPsightsWidget(atoms,electronsVectorCollection);


    singleElectronEnergies(cluster);

    return QApplication::exec();

    /* TODO
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
}
