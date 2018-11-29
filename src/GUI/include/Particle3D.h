//
// Created by heuer on 29.11.18.
//

#ifndef INPSIGHTS_PARTICLE3D_H
#define INPSIGHTS_PARTICLE3D_H

#include <Particle.h>
#include <ElementInfo.h>
#include "Sphere.h"
#include "GuiHelper.h"

template <typename Type>
class Particle3D : public Particle<Type>, public Sphere {
private:
    Q_DISABLE_COPY(Particle3D<Type>)
public:
    Particle3D(Qt3DCore::QEntity *root, Particle<Type> particle)
    : Particle<Type>::Particle(particle),
            Sphere(root,
                   GuiHelper::QColorFromType<Type>(particle.type()),
                   GuiHelper::toQVector3D(particle.position()),
                   GuiHelper::radiusFromType<Type>(particle.type()))
    {
        if(std::is_same<Type,Element>())
            material->setAlpha(0.25f);
        else if (std::is_same<Type,Spin>())
            material->setAlpha(0.5f);
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
