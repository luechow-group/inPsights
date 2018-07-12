//
// Created by Michael Heuer on 20.06.18.
//

#include "Environment.h"
#include "CutoffFunction.h"
#include "ExpansionSettings.h"
#include <SpecialMathFunctions/BoostSphericalHarmonics.h>

SphericalCoordinates::SphericalCoordinates(const Eigen::Vector3d& vec)
        : r(0),theta(0),phi(0) {
    BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, r, theta, phi);
}


Environment::Environment(MolecularGeometry molecularGeometry, Eigen::Vector3d center)
: molecularGeometry_(std::move(molecularGeometry)),
center_(std::move(center)){}

std::vector<std::pair<Particle<int>,SphericalCoordinates>> Environment::selectParticles(int expansionTypeId) const {

    std::vector<std::pair<Particle<int>,SphericalCoordinates>> selectedParticles;

    for (unsigned j = 0; j < unsigned(molecularGeometry_.numberOfEntities()); ++j) {
        const auto &neighbor = molecularGeometry_[j];

        if( neighbor.type() == expansionTypeId || ExpansionSettings::mode == ExpansionSettings::Mode::typeAgnostic) {

            SphericalCoordinates sphericalCoords(neighbor.position()-center_);

            if (CutoffFunction::withinCutoffRadiusQ(sphericalCoords.r)) {
                selectedParticles.emplace_back<std::pair<Particle<int>,SphericalCoordinates>>(
                        {neighbor,sphericalCoords});
            }
        }
    }
    return selectedParticles;
}