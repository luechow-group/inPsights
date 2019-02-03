//
// Created by Michael Heuer on 2019-02-03.
//

#ifndef INPSIGHTS_ENERGYPARTITIONING_H
#define INPSIGHTS_ENERGYPARTITIONING_H

#include <Eigen/Core>

namespace EnergyPartitioning {
    namespace ParticleBased {
        Eigen::VectorXd oneAtomEnergies(const Eigen::MatrixXd &Vnn);

        Eigen::VectorXd oneElectronEnergies(
                const Eigen::VectorXd &Te,
                const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven,
                const Eigen::MatrixXd &Vnn);
    }

}

#endif //INPSIGHTS_ENERGYPARTITIONING_H
