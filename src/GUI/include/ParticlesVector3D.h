//
// Created by Michael Heuer on 02.12.18.
//

#ifndef INPSIGHTS_PARTICLESVECTOR3D_H
#define INPSIGHTS_PARTICLESVECTOR3D_H

#include <Qt3DCore/QEntity>
#include <ParticlesVector.h>
#include "Particle3D.h"

template <typename Type>
class ParticlesVector3D : public ParticlesVector<Type>, public Qt3DCore::QEntity {
public:
    ParticlesVector3D(Qt3DCore::QEntity *root, const ParticlesVector<Type> &particlesVector)
    : ParticlesVector<Type>(particlesVector),
            QEntity(root), particles3D_(0) {

        for (long i = 0; i < ParticlesVector<Type>::numberOfEntities(); ++i) {
            particles3D_.emplace_back(new Particle3D<Type>(this, particlesVector.linkedParticle(i)));
        }
    };

    std::vector<Particle3D<Type>*> particles3D_; // shared pointers?
};

using TypedParticlesVector3D = ParticlesVector3D<int>;
using ElectronsVector3D = ParticlesVector3D<Spin>;
using AtomsVector3D = ParticlesVector3D<Element>;

#endif //INPSIGHTS_PARTICLESVECTOR3D_H
