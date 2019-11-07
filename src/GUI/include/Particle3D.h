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

#ifndef INPSIGHTS_PARTICLE3D_H
#define INPSIGHTS_PARTICLE3D_H

#include <Particle.h>
#include <ElementInfo.h>
#include "Sphere.h"
#include "GuiHelper.h"
#include <memory>

template <typename Type>
class Particle3D : public LinkedParticle<Type>, public Sphere {
public:
    Particle3D(Qt3DCore::QEntity *root, std::shared_ptr<LinkedParticle<Type>> particle)
    : LinkedParticle<Type>::LinkedParticle(*particle),
            Sphere(root,
                   GuiHelper::QColorFromType<Type>(particle->type()),
                   GuiHelper::toQVector3D(particle->position()),
                   GuiHelper::radiusFromType<Type>(particle->type())
                           ) {

        if(std::is_same<Type,Element>())
            material->setAlpha(0.33f);
        else if (std::is_same<Type,Spin>())
            material->setAlpha(0.5f);
    }

    void setPosition(const Eigen::Vector3d &position) override {
        LinkedParticle<Type>::setPosition(position);
        transform->setTranslation(GuiHelper::toQVector3D(position));
    }
};

using TypedParticle3D = Particle3D<int>;
using Electron3D = Particle3D<Spin>;
using Atom3D = Particle3D<Element>;

#endif //INPSIGHTS_PARTICLE3D_H
