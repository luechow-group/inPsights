//
// Created by Michael Heuer on 15.03.18.
//

#include "AtomicNeighborhoodDensityFunction.h"

AtomicNeighborhoodDensityFunction::AtomicNeighborhoodDensityFunction(const PositionsVector& positionsVector, double alpha)
        : alpha_(alpha),
          positionsVector_(positionsVector) {};

double AtomicNeighborhoodDensityFunction::value(const Eigen::Vector3d &rvec) const{

    double value = 0.0;
    for (int i = 0; i < positionsVector_.numberOfEntities(); ++i) {
        double exponent = 0.0;
        exponent += std::pow(rvec.x()-positionsVector_[i].x(),2);
        exponent += std::pow(rvec.y()-positionsVector_[i].y(),2);
        exponent += std::pow(rvec.z()-positionsVector_[i].z(),2);

        value += std::exp(-alpha_*exponent);
    }

    return value;
};
