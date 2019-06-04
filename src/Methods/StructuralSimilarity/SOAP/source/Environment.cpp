//
// Created by Michael Heuer on 20.06.18.
//

#include "Environment.h"
#include "Cutoff.h"
#include "SOAPSettings.h"
#include <SpecialMathFunctions/BoostSphericalHarmonics.h>

using namespace SOAP;

SphericalCoordinates::SphericalCoordinates(const Eigen::Vector3d& vec)
        : r(0),theta(0),phi(0) {
    BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, r, theta, phi);
}

Environment::Environment(MolecularGeometry molecularGeometry, EnumeratedType<int> enumeratedType)
: molecularGeometry_(std::move(molecularGeometry)),
  enumeratedType_(std::move(enumeratedType)){}

std::vector<std::pair<Particle<int>,SphericalCoordinates>> Environment::selectParticles(int expansionTypeId) const {
    auto [foundQ, idx] = molecularGeometry_.findIndexByEnumeratedType(enumeratedType_);
    assert(foundQ && "The enumerated type must exist.");

    auto center = molecularGeometry_[idx].position();

    auto mode = General::settings.mode();
    std::vector<std::pair<Particle<int>,SphericalCoordinates>> selectedParticles;

    for (unsigned i = 0; i < unsigned(molecularGeometry_.numberOfEntities()); ++i) {
        const auto &neighbor = molecularGeometry_[i];

        if( i != idx && (neighbor.type() == expansionTypeId || mode == General::Mode::typeAgnostic) ) {

            SphericalCoordinates sphericalCoords(neighbor.position()-center);

            if (Cutoff::withinCutoffRadiusQ(sphericalCoords.r)) {
                selectedParticles.emplace_back<std::pair<Particle<int>,SphericalCoordinates>>(
                        {neighbor,sphericalCoords});
            }
        }
    }
    return selectedParticles;
}