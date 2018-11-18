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

#include <InPsightsSplashScreen.h>
#include <QHBoxLayout>
#include <QRadioButton>

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
    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);

    setlocale(LC_NUMERIC,"C");
    QWidget* window = new QWidget();

    auto moleculeWidget = new MoleculeWidget();
    auto layout = new QHBoxLayout();
    auto radioButton = new QRadioButton();

    layout->addWidget(moleculeWidget->getWidget(),Qt::AlignLeft);
    layout->addWidget(radioButton);
    window->setLayout(layout);
    window->resize(800,800);

    auto splash = InPsightsSplashScreen::getInPsightsSplashScreen();

    QTimer::singleShot(2000, splash, SLOT(close()));
    QTimer::singleShot(2000, window, SLOT(show()));

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

    singleElectronEnergies(cluster);

    auto electronsVectorCollection = cluster["Structures"].as<std::vector<ElectronsVector>>();
    auto spinCorrelations = cluster["SpinCorrelations"];
    AtomsVector3D(moleculeWidget->getRoot(), atoms);
    ElectronsVector3D(moleculeWidget->getRoot(), electronsVectorCollection[0]);
    //for (const auto & i : electronsVectorCollection)
    //    ElectronsVector3D(root, i);


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
}
