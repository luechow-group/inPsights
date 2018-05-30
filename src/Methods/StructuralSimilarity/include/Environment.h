//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_ENVIRONMENT_H
#define AMOLQCPP_ENVIRONMENT_H

#include <Eigen/Core>
#include <utility>
#include <MolecularGeometry.h>
#include <Type.h>
#include "CutoffFunction.h"
#include <BoostSphericalHarmonics.h>

class SphericalCoordinates{
public:
    SphericalCoordinates(const Eigen::Vector3d& vec)
            : r(0),theta(0),phi(0) {
        BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, r, theta, phi);
    }

    double r, theta, phi;
};

class Environment{
public:
    Environment(MolecularGeometry molecularGeometry, Eigen::Vector3d center)
            : molecularGeometry_(std::move(molecularGeometry)),
              center_(std::move(center)){}

    std::vector<std::pair<Particle<int>,SphericalCoordinates>>
    selectParticles(int expansionTypeId = 0) const {

        std::vector<std::pair<Particle<int>,SphericalCoordinates>> selecetedParticles;

        for (unsigned j = 0; j < unsigned(molecularGeometry_.numberOfEntities()); ++j) {
            const auto &neighbor = molecularGeometry_[j];

            if( neighbor.type() == expansionTypeId || ExpansionSettings::mode == ExpansionMode::Generic) {

                SphericalCoordinates sphericalCoords(neighbor.position()-center_);

                if (CutoffFunction::withinCutoffRadiusQ(sphericalCoords.r)) {
                    selecetedParticles.emplace_back<std::pair<Particle<int>,SphericalCoordinates>>(
                            {neighbor,sphericalCoords});
                }
            }
        }
        return selecetedParticles;
    }

    MolecularGeometry molecularGeometry_;
    Eigen::Vector3d center_; //TODO be careful with many particles located at the same center => getCenterWeight?
};

#endif //AMOLQCPP_ENVIRONMENT_H
