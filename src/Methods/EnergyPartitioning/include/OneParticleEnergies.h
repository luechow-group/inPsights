//
// Created by Michael Heuer on 22.11.18.
//

#ifndef INPSIGHTS_ONEPARTICLEENERGIES_H
#define INPSIGHTS_ONEPARTICLEENERGIES_H

#include <Eigen/Core>
#include <Statistics.h>

class ClusterData;

namespace OneParticleEnergies {

    Eigen::VectorXd oneAtomEnergies(const Eigen::MatrixXd &Vnn);
    Eigen::VectorXd oneElectronEnergies(
            const Eigen::VectorXd &Te,
            const Eigen::MatrixXd &Vee,
            const Eigen::MatrixXd &Ven,
            const Eigen::MatrixXd &Vnn);

    Eigen::VectorXd oneAtomEnergies(const IntraParticlesStatistics& Vnn, const ClusterData &clusterData);
    Eigen::VectorXd oneAtomEnergiesErrors(const IntraParticlesStatistics& Vnn, const ClusterData &clusterData);


    Eigen::VectorXd oneElectronEnergies(const ClusterData &clusterData);
    Eigen::VectorXd oneElectronEnergiesErrors(const ClusterData &clusterData);

}

#endif //INPSIGHTS_ONEPARTICLEENERGIES_H
