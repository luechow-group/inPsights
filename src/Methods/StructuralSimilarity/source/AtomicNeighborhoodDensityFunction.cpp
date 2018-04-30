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
        value += std::exp(-alpha_*(rvec-positionsVector_[i]).squaredNorm());
    }

    return value;
};
