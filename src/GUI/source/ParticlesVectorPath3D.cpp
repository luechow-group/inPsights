/* Copyright (C) 2017-2019 Michael Heuer.
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
