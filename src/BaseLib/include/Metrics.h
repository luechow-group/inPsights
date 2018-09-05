//
// Created by Michael Heuer on 25.04.18.
//

#ifndef AMOLQCPP_METRICS_H
#define AMOLQCPP_METRICS_H


#include <Eigen/Core>
#include "PositionsVector.h"

namespace Metrics{

    double distance(const Eigen::Vector3d& position1,
                    const Eigen::Vector3d& position2);

    double distance(const PositionsVector& positions1,
                    const PositionsVector& positions2);

    Eigen::VectorXd positionDistancesVector(const PositionsVector &positions1, const PositionsVector &positions2);
    
    Eigen::MatrixXd positionalDistances(const PositionsVector& positions);

    Eigen::MatrixXd positionalDistances(const PositionsVector& positions1, 
                                        const PositionsVector& positions2);
}

#endif //AMOLQCPP_METRICS_H
