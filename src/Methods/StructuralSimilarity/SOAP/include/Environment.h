// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ENVIRONMENT_H
#define INPSIGHTS_ENVIRONMENT_H

#include <Eigen/Core>
#include <utility>
#include <MolecularGeometry.h>
#include "SOAPSettings.h"

namespace SOAP {
    class SphericalCoordinates {
    public:
        SphericalCoordinates(const Eigen::Vector3d &vec);

        double r, theta, phi;
    };

    class Environment {
    public:
        Environment(const MolecularGeometry& molecularGeometry, Eigen::Vector3d  center);
        Environment(const MolecularGeometry& molecularGeometry, EnumeratedType<int> enumeratedType);

        std::vector<std::pair<Particle<int>, SphericalCoordinates>> selectParticles(int expansionTypeId) const;
        bool selectParticleQ(unsigned index, const TypedParticle& neighbor, int expansionTypeId) const;

    private:
        static long getOwnIndex(const MolecularGeometry& molecularGeometry, EnumeratedType<int> enumeratedType) ;

        const MolecularGeometry& molecularGeometry_;
        long ownIdx_;
        bool agnosticQ_;
        Eigen::Vector3d  center_;
    };
}

#endif //INPSIGHTS_ENVIRONMENT_H
