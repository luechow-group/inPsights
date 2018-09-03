//
// Created by Michael Heuer on 25.04.18.
//

#include "Metrics.h"

double Metrics::distance(const Eigen::Vector3d& position1, const Eigen::Vector3d& position2){
    return (position1-position2).norm();
}


double Metrics::distance(const PositionsVector& positions1,
                         const PositionsVector& positions2){
    assert(positions1.numberOfEntities() == positions2.numberOfEntities()
           && "Both PositionVectors must have the same size.");
    return (positions1.asEigenVector()-positions2.asEigenVector()).norm();
}

Eigen::MatrixXd Metrics::positionalDistances(const PositionsVector &positions){
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(positions.numberOfEntities(),positions.numberOfEntities());
    for (size_t i = 0; i < d.rows(); i++)
        for (size_t j = i + 1; j < d.cols(); j++)
            d(i,j) = Metrics::distance(positions[i], positions[j]);

    // symmetrization
    return d.selfadjointView<Eigen::Upper>();
};

Eigen::MatrixXd Metrics::positionalDistances(const PositionsVector &positions1, const PositionsVector &positions2){
    assert(positions1.numberOfEntities() == positions2.numberOfEntities());

    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(positions1.numberOfEntities(),positions1.numberOfEntities());

    for (size_t i = 0; i < d.rows(); i++)
        for (size_t j = 0; j < d.cols(); j++)
            d(i,j) = Metrics::distance(positions1[i], positions2[j]);

    return d;
};
