//
// Created by heuer on 15.11.18.
//

#include <iostream>
#include <ParticlesVector.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <InPsightsWidget.h>
#include <QApplication>
#include <yaml-cpp/yaml.h>

void singleElectronEnergies(const YAML::Node &cluster) {
    auto Te = cluster["Te"];
    auto Vee = cluster["Vee"];
    auto Ven = cluster["Ven"];
    auto nElectrons = Te.size();
    auto nAtoms = Ven.begin()->second.size();
    std::vector<double> Ve(nElectrons);
    std::vector<double> VeErr(nElectrons);

    for (unsigned long i = 0; i < nElectrons; ++i) {
        Ve[i] = Te[i][0].as<double>();
        VeErr[i] = Te[i][1].as<double>();

        for (unsigned long  k = 0; k < nAtoms; ++k){
            Ve[i] += Ven[i][k][0].as<double>();
            VeErr[i] += Ven[i][k][1].as<double>();
        }

        for (unsigned long  j = i+1; j < nElectrons; ++j){
            Ve[i] += 0.5*Vee[i][j][0].as<double>();
            VeErr[i] += 0.5*Vee[i][j][1].as<double>();
        }
        //TODO FEHLERADDITION
        std::cout << "Ee = " << Ve[i] << "+/-" << VeErr[i] << "  (ERROR is wrong -> add error propagation)" << std::endl;
    }
};


int main(int argc, char *argv[]) {

    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    new InPsightsWidget();

    //singleElectronEnergies(cluster);

    return QApplication::exec();

    /* TODO
     * METHOD
     * - spatial permutations (by value range or struct sim)
     *
     * GUI
     *  - dotted and dashed lines
     * - select clusters
     * - show energies
     * - show eigenvectors (.wf needed)
     */
}
