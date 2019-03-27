//
// Created by Michael Heuer on 17.05.18.
//

#ifndef INPSIGHTS_SINKHORN_H
#define INPSIGHTS_SINKHORN_H

#include <Eigen/Core>
#include "ExpansionSettings.h"
namespace Sinkhorn{

    Eigen::MatrixXd Pgamma(const Eigen::MatrixXd &C,
                           double gamma = SOAPExpansion::settings.gamma(),
                           double eps = std::numeric_limits<double>::epsilon());

    double distance(Eigen::MatrixXd correlationMatrix, double gamma = SOAPExpansion::settings.gamma(), double eps = std::numeric_limits<double>::epsilon());
}

#endif //INPSIGHTS_SINKHORN_H
