//
// Created by Michael Heuer on 22.11.18.
//

#include <OneParticleEnergies.h>
#include <ClusterData.h>
#include <math.h>

Eigen::VectorXd
OneParticleEnergies::oneAtomEnergies(const IntraParticlesStatistics& Vnn, const ClusterData &clusterData) {
    //auto nElectrons = clusterData.VenStats_.rows(); //TODO delete?
    auto nAtoms = clusterData.VenStats_.cols();
    //const auto& Ven = clusterData.VenStats_; //TODO delete?

    Eigen::VectorXd En = Eigen::VectorXd::Zero(nAtoms);

    for (Eigen::Index k = 0; k < nAtoms; ++k) {
        for (Eigen::Index l = k + 1; l < nAtoms; ++l) {
            En[k] += 0.5 * Vnn.mean()(k, l);
            En[l] += 0.5 * Vnn.mean()(k, l);
        }

        //for (Eigen::Index  i = 0; i < nElectrons; ++i)
        //    En[k] += 0.5 * Ven.mean()(i,k); //TODO delete?
    }
    return En;
}

Eigen::VectorXd
OneParticleEnergies::oneAtomEnergiesErrors(const IntraParticlesStatistics &Vnn, const ClusterData &clusterData) {
    //auto nElectrons = clusterData.VenStats_.rows();
    auto nAtoms = clusterData.VenStats_.cols();
    //const auto& Ven = clusterData.VenStats_;

    Eigen::VectorXd EnErr = Eigen::VectorXd::Zero(nAtoms);

    for (Eigen::Index k = 0; k < nAtoms; ++k) {
        for (Eigen::Index l = k + 1; l < nAtoms; ++l) {
            EnErr[k] = 0.5 * std::sqrt(std::pow(EnErr[k], 2) + std::pow(Vnn.standardError()(k, l), 2));
            EnErr[l] = 0.5 * std::sqrt(std::pow(EnErr[l], 2) + std::pow(Vnn.standardError()(k, l), 2));
        }

        //for (Eigen::Index  i = 0; i < nElectrons; ++i)
        //    EnErr[k] = 0.5 * std::sqrt(std::pow(EnErr[k], 2) + std::pow(Ven.standardError()(i,k), 2)); //TODO delete?
    }
    return EnErr;
}

Eigen::VectorXd OneParticleEnergies::oneElectronEnergies(const ClusterData &clusterData) {
    auto nElectrons = clusterData.VenStats_.rows();
    auto nAtoms = clusterData.VenStats_.cols();
    const auto& Te = clusterData.TeStats_;
    const auto& Vee = clusterData.VeeStats_;
    const auto& Ven = clusterData.VenStats_;

    Eigen::VectorXd Ee = Eigen::VectorXd::Zero(nElectrons);

    for (Eigen::Index i = 0; i < nElectrons; ++i) {
        Ee[i] += Te.mean()[i];

        for (Eigen::Index k = 0; k < nAtoms; ++k)
            Ee[i] += Ven.mean()(i, k);
            //Ee[i] += 0.5 * Ven.mean()(i, k); //TODO delete?

        for (Eigen::Index  j = i + 1; j < nElectrons; ++j) {
            Ee[i] += 0.5 * Vee.mean()(i, j);
            Ee[j] += 0.5 * Vee.mean()(i, j);
        }
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
        EeErr[i] += std::pow(Te.standardError()[i], 2);

        for (Eigen::Index  k = 0; k < nAtoms; ++k)
            EeErr[i] +=  std::pow(Ven.standardError()(i,k), 2);
            //EeErr[i] = 0.5*std::sqrt(std::pow(EeErr[i], 2) + std::pow(Ven.standardError()(i,k), 2)); //TODO delete?

        for (Eigen::Index  j = i + 1; j < nElectrons; ++j) {
            EeErr[i] += std::pow(0.5*Vee.standardError()(i, j), 2);
            EeErr[j] += std::pow(0.5*Vee.standardError()(i, j), 2);
        }
    }

    EeErr.cwiseSqrt();

    return EeErr;
}
