//
// Created by Michael Heuer on 22.11.18.
//

#include <OneParticleEnergies.h>
#include <ClusterData.h>
#include <math.h>

Eigen::VectorXd
OneParticleEnergies::oneAtomEnergies(const IntraParticlesStatistics& Vnn, const ClusterData &clusterData) {
    auto nElectrons = clusterData.VenStats_.rows();
    auto nAtoms = clusterData.VenStats_.cols();
    const auto& Ven = clusterData.VenStats_;

    Eigen::VectorXd En = Eigen::VectorXd::Zero(nAtoms);

    for (Eigen::Index k = 0; k < nAtoms; ++k) {
        for (Eigen::Index l = k + 1; l < nAtoms; ++l)
            En[k] += 0.5 * Vnn.mean()(k,l);

        for (Eigen::Index  i = 0; i < nElectrons; ++i)
            En[k] += 0.5 * Ven.mean()(i,k);
    }
    return En;
}

Eigen::VectorXd
OneParticleEnergies::oneAtomEnergiesErrors(const IntraParticlesStatistics &Vnn, const ClusterData &clusterData) {
    auto nElectrons = clusterData.VenStats_.rows();
    auto nAtoms = clusterData.VenStats_.cols();
    const auto& Ven = clusterData.VenStats_;

    Eigen::VectorXd EnErr = Eigen::VectorXd::Zero(nAtoms);

    for (Eigen::Index i = 0; i < nAtoms; ++i) {
        for (Eigen::Index j = i + 1; j < nAtoms; ++j)
            EnErr[i] = 0.5 * std::sqrt(std::pow(EnErr[i], 2) + std::pow(Vnn.standardError()(i,j), 2));

        for (Eigen::Index  k = 0; k < nElectrons; ++k)
            EnErr[i] = 0.5 * std::sqrt(std::pow(EnErr[i], 2) + std::pow(Ven.standardError()(k,i), 2));
    }
    return EnErr;
}

Eigen::VectorXd OneParticleEnergies::oneElectronEnergies(const ClusterData &clusterData) {
    auto nElectrons = clusterData.VenStats_.rows();
    auto nAtoms = clusterData.VenStats_.cols();
    const auto& Te = clusterData.TeStats_;
    const auto& Vee = clusterData.VeeStats_;
    const auto& Ven = clusterData.VenStats_;

    Eigen::VectorXd Ee= Eigen::VectorXd::Zero(nElectrons);

    for (Eigen::Index i = 0; i < nElectrons; ++i) {
        Ee[i] = Te.mean()[i];

        for (Eigen::Index  k = 0; k < nAtoms; ++k)
            Ee[i] += 0.5*Ven.mean()(i,k);

        for (Eigen::Index  j = i + 1; j < nElectrons; ++j)
            Ee[i] += 0.5 * Vee.mean()(i,j);
    }
    return Ee;
}

Eigen::VectorXd OneParticleEnergies::oneElectronEnergiesErrors(const ClusterData &clusterData) {
    auto nElectrons = clusterData.VenStats_.rows();
    auto nAtoms = clusterData.VenStats_.cols();
    const auto& Te = clusterData.TeStats_;
    const auto& Vee = clusterData.VeeStats_;
    const auto& Ven = clusterData.VenStats_;

    Eigen::VectorXd EeErr = Eigen::VectorXd::Zero(nElectrons);

    for (Eigen::Index i = 0; i < nElectrons; ++i) {
        EeErr[i] = Te.standardError()[i];

        for (Eigen::Index  k = 0; k < nAtoms; ++k)
            EeErr[i] = 0.5*std::sqrt(std::pow(EeErr[i], 2) + std::pow(Ven.standardError()(i,k), 2));

        for (Eigen::Index  j = i + 1; j < nElectrons; ++j)
            EeErr[i] = 0.5 * std::sqrt(std::pow(EeErr[i], 2) + std::pow(Vee.standardError()(i,j), 2));
    }
    return EeErr;
}
