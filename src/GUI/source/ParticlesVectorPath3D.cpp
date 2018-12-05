//
// Created by Michael Heuer on 12.11.17.
//

#include <Polyline.h>
#include <ParticlesVectorPath3D.h>
#include <Particle3D.h>
#include <GuiHelper.h>

ParticlesVectorPath3D::ParticlesVectorPath3D(Qt3DCore::QEntity *root,
                                                   const ElectronsVectorCollection &electronsVectorCollection,
                                                   float radius) {


    auto numberOfParticles = electronsVectorCollection[0].numberOfEntities();
    std::vector<std::vector<QVector3D>> pointsList(numberOfParticles);
    for (long i = 0; i < numberOfParticles; ++i) { // iterate over particles

        for (long j = 0; j < electronsVectorCollection.numberOfEntities(); ++j)
            pointsList[i].emplace_back(GuiHelper::toQVector3D(electronsVectorCollection[j][i].position()));

        auto spinType = electronsVectorCollection.typesVector()[i];

        if (spinType == Spin::alpha) {
            new Polyline(root,GuiHelper::QColorFromType<Spin>(Spin::alpha) , pointsList[i], radius);
        }
        else {
            new Polyline(root,GuiHelper::QColorFromType<Spin>(Spin::beta) , pointsList[i], radius);
        }

    }

}
