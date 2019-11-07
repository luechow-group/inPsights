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
            QEntity(root),
            connections_(new Qt3DCore::QEntity(this)),
            correlations_(new Qt3DCore::QEntity(this)),
            particles3D_(0) {

        drawParticles(true);
    };

public slots:
    void drawParticles(bool drawQ = true) {
        if(drawQ && particles3D_.empty())
            for(long i = 0; i < ParticlesVector<Type>::numberOfEntities(); ++i)
                particles3D_.emplace_back(new Particle3D<Type>(this, this->linkedParticle(i)));
        else
            for(auto it = particles3D_.begin(); it != particles3D_.end(); it++)
                particles3D_.erase(it);
};

    void drawConnections() {};

    void deleteConnections(){
        connections_->deleteLater();
        connections_ = new Qt3DCore::QEntity(this);
    }

    void deleteCorrelations() {
        correlations_->deleteLater();
        correlations_ = new Qt3DCore::QEntity(this);
    }

    Qt3DCore::QEntity* connections_, *correlations_;
    std::vector<std::shared_ptr<Particle3D<Type>>> particles3D_;
};

using AtomsVector3D = ParticlesVector3D<Element>;
using ElectronsVector3D = ParticlesVector3D<Spin>;

template<>
void AtomsVector3D::drawConnections();
template<>
void ElectronsVector3D::drawConnections();


#endif //INPSIGHTS_PARTICLESVECTOR3D_H
