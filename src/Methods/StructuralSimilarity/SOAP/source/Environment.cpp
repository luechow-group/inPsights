// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Environment.h"
#include "Cutoff.h"
#include <SpecialMathFunctions/BoostSphericalHarmonics.h>

#include <utility>

using namespace SOAP;

SphericalCoordinates::SphericalCoordinates(const Eigen::Vector3d &vec)
        : r(0), theta(0), phi(0) {
    BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, r, theta, phi);
}

Environment::Environment(const MolecularGeometry &molecularGeometry, Eigen::Vector3d center)
        : molecularGeometry_(molecularGeometry),
          ownIdx_(-1), // this influences the behaviour of selectParticle()
          agnosticQ_(General::settings.mode() == General::Mode::typeAgnostic),
          center_(std::move(center)) {}

Environment::Environment(const MolecularGeometry &molecularGeometry, EnumeratedType<int> enumeratedType)
        : molecularGeometry_(molecularGeometry),
          ownIdx_(getOwnIndex(molecularGeometry, enumeratedType)),
          agnosticQ_(General::settings.mode() == General::Mode::typeAgnostic),
          center_(molecularGeometry.positions()[ownIdx_]) {}

long Environment::getOwnIndex(const MolecularGeometry &molecularGeometry, EnumeratedType<int> enumeratedType) {
    assert(molecularGeometry.findIndexByEnumeratedType(enumeratedType).first && "The enumerated type must exist.");
    return molecularGeometry.findIndexByEnumeratedType(enumeratedType).second;
}

std::vector<std::pair<Particle<int>, SphericalCoordinates>> Environment::selectParticles(int expansionTypeId) const {
    std::vector<std::pair<Particle<int>, SphericalCoordinates>> selectedParticles;

    for (unsigned i = 0; i < unsigned(molecularGeometry_.numberOfEntities()); ++i) {
        const auto &neighbor = molecularGeometry_[i];

        if (selectParticleQ(i, neighbor, expansionTypeId)) {
            SphericalCoordinates sphericalCoords(neighbor.position() - center_);

            if (Cutoff::withinCutoffRadiusQ(sphericalCoords.r)) {
                selectedParticles.emplace_back<std::pair<Particle<int>, SphericalCoordinates>>(
                        {neighbor, sphericalCoords});
            }
        }
    }
    return selectedParticles;
}

bool Environment::selectParticleQ(unsigned int index, const TypedParticle &neighbor, int expansionTypeId) const {
    if(neighbor.type() == expansionTypeId || agnosticQ_) // type is ok
        return index != ownIdx_; // not the same particle // TODO change this?
    else
        return false;
}
