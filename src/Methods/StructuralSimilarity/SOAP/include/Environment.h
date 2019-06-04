//
// Created by Michael Heuer on 03.05.18.
//

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
