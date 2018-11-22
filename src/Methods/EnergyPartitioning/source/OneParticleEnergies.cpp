//
// Created by Michael Heuer on 22.11.18.
//

#include <OneParticleEnergies.h>

void OneParticleEnergies::oneAtomEnergies(const YAML::Node &cluster, const YAML::Node &Vnn) {

    auto Ven = cluster["Ven"];
    auto nElectrons = Ven.size();
    auto nAtoms = Ven.begin()->second.size();
    std::vector<double> Vn(nAtoms);
    std::vector<double> VnErr(nAtoms);


    for (unsigned long i = 0; i < nAtoms; ++i) {

        for (unsigned long j = i + 1; j < nAtoms; ++j) {
            Vn[i] += 0.5 * Vnn[i][j][0].as<double>();
            VnErr[i] = 0.5 * std::sqrt(std::pow(VnErr[i], 2) + std::pow(Vnn[i][j][1].as<double>(), 2));
        }

        for (unsigned long k = 0; k < nElectrons; ++k) {
            Vn[i] += 0.5 * Ven[k][i][0].as<double>();
            VnErr[i] = 0.5 * std::sqrt(std::pow(VnErr[i], 2) + std::pow(Ven[k][i][1].as<double>(), 2));
        }

        std::cout << "En = " << Vn[i] << "+/-" << VnErr[i] << std::endl;
    }
}

void OneParticleEnergies::oneElectronEnergies(const YAML::Node &cluster) {
    auto Te = cluster["Te"];
    auto Vee = cluster["Vee"];
    auto Ven = cluster["Ven"];
    auto nElectrons = Ven.size();
    auto nAtoms = Ven.begin()->second.size();
    std::vector<double> Ve(nElectrons);
    std::vector<double> VeErr(nElectrons);

    for (unsigned long i = 0; i < nElectrons; ++i) {
        Ve[i] = Te[i][0].as<double>();
        VeErr[i] = Te[i][1].as<double>();

        for (unsigned long k = 0; k < nAtoms; ++k) {
            Ve[i] += 0.5*Ven[i][k][0].as<double>();
            VeErr[i] = 0.5*std::sqrt(std::pow(VeErr[i], 2) + std::pow(Ven[i][k][1].as<double>(), 2));
        }

        for (unsigned long j = i + 1; j < nElectrons; ++j) {
            Ve[i] += 0.5 * Vee[i][j][0].as<double>();
            VeErr[i] = 0.5 * std::sqrt(std::pow(VeErr[i], 2) + std::pow(Vee[i][j][1].as<double>(), 2));
        }
        std::cout << "Ee = " << Ve[i] << "+/-" << VeErr[i] << std::endl;
    }
};
