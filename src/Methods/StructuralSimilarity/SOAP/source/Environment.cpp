/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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