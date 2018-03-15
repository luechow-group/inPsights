//
// Created by Michael Heuer on 15.03.18.
//

#ifndef AMOLQCPP_ATOMICNEIGHBORHOODDENSITYFUNCTION_H
#define AMOLQCPP_ATOMICNEIGHBORHOODDENSITYFUNCTION_H

#include "LebedevSphericalIntegration/SpatialFunction.h"
#include "PositionsVector.h"

class AtomicNeighborhoodDensityFunction : public SpatialFunction{

public:
    explicit AtomicNeighborhoodDensityFunction(const PositionsVector& positionsVector, double alpha = 1.0);

    double value(const Eigen::Vector3d &rvec) const override;

private:
    PositionsVector positionsVector_;
    double alpha_;
};

#endif //AMOLQCPP_ATOMICNEIGHBORHOODDENSITYFUNCTION_H
