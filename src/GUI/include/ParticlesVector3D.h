// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PARTICLESVECTOR3D_H
#define INPSIGHTS_PARTICLESVECTOR3D_H

#include <Qt3DCore/QEntity>
#include <ParticlesVector.h>
#include "Particle3D.h"
#include "ColorPalette.h"
#include <Metrics.h>

template <typename Type>
class ParticlesVector3D : public ParticlesVector<Type>, public Qt3DCore::QEntity {
public:
    ParticlesVector3D(Qt3DCore::QEntity *root, const ParticlesVector<Type> &particlesVector, bool coloredQ = false)
    : ParticlesVector<Type>(particlesVector),
            QEntity(root),
            connections_(new Qt3DCore::QEntity(this)),
            correlations_(new Qt3DCore::QEntity(this)),
            coloredQ_(coloredQ),
            particles3D_(0) {
        drawParticles(true);
    };

public slots:
    void drawParticles(bool drawQ = true) {
        if(drawQ && particles3D_.empty())
            for(long i = 0; i < ParticlesVector<Type>::numberOfEntities(); ++i)
                if(coloredQ_) {
                    particles3D_.emplace_back(new Particle3D<Type>(this, this->operator[](i),
                                              ColorPalette::colorFunction(i), 1.0f/3.0f, 2, 4, 1.0));
                } else {
                    bool violet = false;
                    for(long j = 0; j < ParticlesVector<Type>::numberOfEntities(); ++j)
                        if (j != i){
                            if (Metrics::distance(this->operator[](i).position(), this->operator[](j).position()) < 0.01){
                                violet = true;
                                break;
                            }
                        }
                    if (not violet) {
                        particles3D_.emplace_back(new Particle3D<Type>(this, this->operator[](i)));
                    }
                    else {
                        particles3D_.emplace_back(new Particle3D<Type>(this, this->operator[](i), Qt::darkMagenta));
                    }
                }
        else
            for(auto it = particles3D_.begin(); it != particles3D_.end(); it++)
                particles3D_.erase(it);
};

    void drawConnections(const double &limit) {};

    void deleteConnections(){
        connections_->deleteLater();
        connections_ = new Qt3DCore::QEntity(this);
    }

    void deleteCorrelations() {
        correlations_->deleteLater();
        correlations_ = new Qt3DCore::QEntity(this);
    }

    Qt3DCore::QEntity* connections_, *correlations_;
    bool coloredQ_;
    std::vector<std::shared_ptr<Particle3D<Type>>> particles3D_;
};

using AtomsVector3D = ParticlesVector3D<Element>;
using ElectronsVector3D = ParticlesVector3D<Spin>;

template<>
void AtomsVector3D::drawConnections(const double &limit);
template<>
void ElectronsVector3D::drawConnections(const double &limit);


#endif //INPSIGHTS_PARTICLESVECTOR3D_H
