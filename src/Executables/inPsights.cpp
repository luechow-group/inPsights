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
#include <QGroupBox>
#include <QSpinBox>
#include <QLabel>

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

    QWidget* inPsightsWidget = new QWidget();

    auto moleculeWidget = new MoleculeWidget();
    auto layout = new QHBoxLayout();
    QSpinBox *integerSpinBox = new QSpinBox;

    layout->addWidget(moleculeWidget, Qt::AlignLeft);
    layout->addWidget(integerSpinBox);

    inPsightsWidget->setLayout(layout);
    inPsightsWidget->resize(800,800);

    auto splash = InPsightsSplashScreen::getInPsightsSplashScreen();

    YAML::Node doc = YAML::LoadFile("raw.yml");

    auto nClusters =doc["Clusters"].size();

    auto spinBoxesGroup = new QGroupBox();
    QLabel *integerLabel = new QLabel("Cluster Id:");
    integerSpinBox->setRange(0, nClusters-1);
    integerSpinBox->setSingleStep(1);
    integerSpinBox->setValue(0);

    size_t clusterId = 0;
    auto cluster = doc["Clusters"][clusterId];


    QTimer::singleShot(2000, splash, SLOT(close()));
    QTimer::singleShot(2000, inPsightsWidget, SLOT(show()));
    moleculeWidget->show();

    std::string resultFilename;

    if (resultFilename.empty()) {
        bool inputArgumentsFoundQ = handleCommandlineArguments(argc, argv, resultFilename, clusterId);
        if (!inputArgumentsFoundQ) return 1;
    }



    singleElectronEnergies(cluster);

    auto electronsVectorCollection = cluster["Structures"].as<std::vector<ElectronsVector>>();
    auto spinCorrelations = cluster["SpinCorrelations"];

    AtomsVector3D(moleculeWidget->getRoot(), doc["Atoms"].as<AtomsVector>());
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
