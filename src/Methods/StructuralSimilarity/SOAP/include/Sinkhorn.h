//
// Created by Michael Heuer on 17.05.18.
//

#ifndef AMOLQCPP_SINKHORN_H
#define AMOLQCPP_SINKHORN_H

#include <Eigen/Core>
#include "ExpansionSettings.h"
namespace Sinkhorn{

    Eigen::MatrixXd Pgamma(const Eigen::MatrixXd &C,
                           double gamma = ExpansionSettings::gamma,
                           double eps = std::numeric_limits<double>::epsilon());

    double distance(Eigen::MatrixXd correlationMatrix, double gamma = ExpansionSettings::gamma, double eps = std::numeric_limits<double>::epsilon());
}

#endif //AMOLQCPP_SINKHORN_H