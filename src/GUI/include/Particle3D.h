//
// Created by heuer on 29.11.18.
//

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
