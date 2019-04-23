//
// Created by Michael Heuer on 18.04.18.
//

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
