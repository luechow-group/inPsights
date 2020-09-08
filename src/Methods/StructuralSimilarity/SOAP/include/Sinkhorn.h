// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SINKHORN_H
#define INPSIGHTS_SINKHORN_H

#include <Eigen/Core>
#include "SOAPSettings.h"


namespace SOAP {
    namespace Sinkhorn {
        Eigen::MatrixXd Pgamma(const Eigen::MatrixXd &C,
                               double gamma = General::settings.sinkhornGamma(),
                               double eps = std::numeric_limits<double>::epsilon());

        double distance(const Eigen::MatrixXd& correlationMatrix, double gamma = General::settings.sinkhornGamma(),
                        double eps = std::numeric_limits<double>::epsilon());
    }
}

#endif //INPSIGHTS_SINKHORN_H
