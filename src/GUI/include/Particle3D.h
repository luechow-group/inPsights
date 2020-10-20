// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PARTICLE3D_H
#define INPSIGHTS_PARTICLE3D_H

#include <Particle.h>
#include <ElementInfo.h>
#include "Sphere.h"
#include "GuiHelper.h"
#include <memory>

template <typename Type>
class Particle3D : public Particle<Type>, public Sphere {
public:
    Particle3D(Qt3DCore::QEntity *root, const Particle<Type>& particle)
    : Particle<Type>::Particle(particle),
            Sphere(root,
                   GuiHelper::QColorFromType<Type>(this->type()),
                   GuiHelper::toQVector3D(this->position()),
                   GuiHelper::radiusFromType<Type>(this->type()),
                   1.0, 12, 24
                           ) {
        if(std::is_same<Type,Element>())
            material->setAlpha(0.33f);
        else if (std::is_same<Type,Spin>())
            material->setAlpha(0.5f);
    }

    // colored constructor
    Particle3D(Qt3DCore::QEntity *root, const Particle<Type>& particle, const QColor& color)
            : Particle<Type>::Particle(particle),
              Sphere(root, color,
                     GuiHelper::toQVector3D(this->position()),
                     GuiHelper::radiusFromType<Type>(this->type())/3.0,
                     1.0, 2, 4
              ) {
        if(std::is_same<Type,Element>())
            material->setAlpha(0.33);
        else if (std::is_same<Type,Spin>())
            material->setAlpha(1.0f);
    }


    void setPosition(const Eigen::Vector3d &position) override {
        Particle<Type>::setPosition(position);
        transform->setTranslation(GuiHelper::toQVector3D(position));
    }
};


using TypedParticle3D = Particle3D<int>;
using Electron3D = Particle3D<Spin>;
using Atom3D = Particle3D<Element>;

#endif //INPSIGHTS_PARTICLE3D_H
