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

#ifndef INPSIGHTS_ENVIRONMENT_H
#define INPSIGHTS_ENVIRONMENT_H

#include <Eigen/Core>
#include <utility>
#include <MolecularGeometry.h>

namespace SOAP {
    class SphericalCoordinates {
    public:
        SphericalCoordinates(const Eigen::Vector3d &vec);

        double r, theta, phi;
    };

    class Environment {
    public:
        Environment(MolecularGeometry molecularGeometry, EnumeratedType<int> enumeratedType);

        std::vector<std::pair<Particle<int>, SphericalCoordinates>> selectParticles(int expansionTypeId = 0) const;

        MolecularGeometry molecularGeometry_;
        EnumeratedType<int> enumeratedType_;
    };
}

#endif //INPSIGHTS_ENVIRONMENT_H
