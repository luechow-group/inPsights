//
// Created by Michael Heuer on 22.11.18.
//

#include <OneParticleEnergies.h>

Eigen::VectorXd OneParticleEnergies::oneAtomEnergies(const Eigen::MatrixXd &Vnn) {
    auto nAtoms = Vnn.rows();

    Eigen::VectorXd En = Eigen::VectorXd::Zero(nAtoms);

    for (Eigen::Index k = 0; k < nAtoms; ++k) {
        for (Eigen::Index l = k + 1; l < nAtoms; ++l) {
            En[k] += 0.5 * Vnn(k, l);
            En[l] += 0.5 * Vnn(k, l);
        }
    }
    return En;
}

Eigen::VectorXd OneParticleEnergies::oneElectronEnergies(const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                         const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
    auto nElectrons = Ven.rows();
    auto nAtoms = Ven.cols();

    Eigen::VectorXd Ee = Eigen::VectorXd::Zero(nElectrons);

    for (Eigen::Index i = 0; i < nElectrons; ++i) {
        Ee[i] += Te[i];

        for (Eigen::Index k = 0; k < nAtoms; ++k)
            Ee[i] += Ven(i, k);

        for (Eigen::Index  j = i + 1; j < nElectrons; ++j) {
            Ee[i] += 0.5 * Vee(i, j);
            Ee[j] += 0.5 * Vee(i, j);
        }
    }
    return Ee;
}

