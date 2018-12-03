//
// Created by Michael Heuer on 02.12.18.
//

#ifndef INPSIGHTS_PARTICLESVECTOR3D_H
#define INPSIGHTS_PARTICLESVECTOR3D_H

#include <Qt3DCore/QEntity>
#include <ParticlesVector.h>
#include "Particle3D.h"
#include <IConnection.h>

template <typename Type>
class ParticlesVector3D : public ParticlesVector<Type>, public Qt3DCore::QEntity {
public:
    ParticlesVector3D(Qt3DCore::QEntity *root, const ParticlesVector<Type> &particlesVector)
    : ParticlesVector<Type>(particlesVector),
            QEntity(root),
            connections_(new Qt3DCore::QEntity(this)),
            correlations_(new Qt3DCore::QEntity(this)),
            particles3D_(0) {

        for (long i = 0; i < ParticlesVector<Type>::numberOfEntities(); ++i) {
            particles3D_.emplace_back(new Particle3D<Type>(this, particlesVector.linkedParticle(i)));
        }
    };

    void drawConnections() {};

    void deleteConnections(){
        connections_->deleteLater();
        connections_ = new Qt3DCore::QEntity(this);
    }

    void deleteCorrelations() {
        std::cout << "delete corr" << std::endl;
        correlations_->deleteLater();
        correlations_ = new Qt3DCore::QEntity(this);
    }

    Qt3DCore::QEntity* connections_, *correlations_;
    std::vector<std::shared_ptr<Particle3D<Type>>> particles3D_;
    std::vector<IConnection*> iConnections_;
};

using AtomsVector3D = ParticlesVector3D<Element>;
using ElectronsVector3D = ParticlesVector3D<Spin>;

template<>
void AtomsVector3D::drawConnections();
template<>
void ElectronsVector3D::drawConnections();


#endif //INPSIGHTS_PARTICLESVECTOR3D_H
