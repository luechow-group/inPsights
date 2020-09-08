// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ANGULARBASIS_H
#define INPSIGHTS_ANGULARBASIS_H

#include <Eigen/Core>

namespace SOAP {
    namespace AngularBasis {
        std::complex<double> computeCoefficient(unsigned l, int m, const Eigen::Vector3d &position);

        std::complex<double> computeCoefficient(unsigned l, int m, double theta, double phi);
    }
}
#endif //INPSIGHTS_ANGULARBASIS_H
