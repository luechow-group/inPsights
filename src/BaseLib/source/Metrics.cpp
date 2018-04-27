//
// Created by Michael Heuer on 25.04.18.
//


#include "Metrics.h"

double Metrics::distance(const Eigen::Vector3d& position1, const Eigen::Vector3d& position2){
    return (position1-position2).norm();
}


double Metrics::distance(const PositionsVector& positionsVector1,
                         const PositionsVector& positionsVector2){
    assert(positionsVector1.numberOfEntities() == positionsVector2.numberOfEntities()
           && "Both PositionVectors must have the same size.");
    return (positionsVector1.positionsAsEigenVector()-positionsVector2.positionsAsEigenVector()).norm();
}