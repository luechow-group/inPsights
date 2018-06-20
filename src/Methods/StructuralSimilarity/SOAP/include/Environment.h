//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_ENVIRONMENT_H
#define AMOLQCPP_ENVIRONMENT_H

#include <Eigen/Core>
#include <utility>
#include <MolecularGeometry.h>

class SphericalCoordinates{
public:
    SphericalCoordinates(const Eigen::Vector3d& vec);
    double r, theta, phi;
};

class Environment{
public:
    Environment(MolecularGeometry molecularGeometry, Eigen::Vector3d center);

    std::vector<std::pair<Particle<int>,SphericalCoordinates>> selectParticles(int expansionTypeId = 0) const;

    MolecularGeometry molecularGeometry_;
    Eigen::Vector3d center_; //TODO be careful with many particles located at the same center => getCenterWeight?
};

#endif //AMOLQCPP_ENVIRONMENT_H
