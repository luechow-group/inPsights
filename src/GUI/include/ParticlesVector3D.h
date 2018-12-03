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
            QEntity(root), connections_(this), particles3D_(0) {

        for (long i = 0; i < ParticlesVector<Type>::numberOfEntities(); ++i) {
            particles3D_.emplace_back(new Particle3D<Type>(this, particlesVector.linkedParticle(i)));
        }
    };

    ~ParticlesVector3D() { /*QT manages destruction*/ };

    // TODO add to vecotr of IConnections
    void drawConnections(){
        //dummy
    }

    void deleteConnections(){
        connections_->deleteLater();
    }

    Qt3DCore::QEntity* connections_;
    std::vector<Particle3D<Type>*> particles3D_; // shared pointers?
};

using AtomsVector3D = ParticlesVector3D<Element>;
using ElectronsVector3D = ParticlesVector3D<Spin>;

template<>
void AtomsVector3D::drawConnections();
template<>
void ElectronsVector3D::drawConnections();


#endif //INPSIGHTS_PARTICLESVECTOR3D_H
