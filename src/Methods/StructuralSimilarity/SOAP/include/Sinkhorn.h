/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_SINKHORN_H
#define INPSIGHTS_SINKHORN_H

#include <Eigen/Core>
#include "SOAPSettings.h"


namespace SOAP {
    namespace Sinkhorn {

        Eigen::MatrixXd Pgamma(const Eigen::MatrixXd &C,
                               double gamma = General::settings.sinkhornGamma(),
                               double eps = std::numeric_limits<double>::epsilon());

        double distance(Eigen::MatrixXd correlationMatrix, double gamma = General::settings.sinkhornGamma(),
                        double eps = std::numeric_limits<double>::epsilon());
    }
}

#endif //INPSIGHTS_SINKHORN_H
