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

#include <Bonds3D.h>
#include <Bond3D.h>
#include <Metrics.h>

Bonds3D::Bonds3D(AtomsVector3D *atomsVector3D, double bondDrawingLimit )
        : IConnection(atomsVector3D->connections_),
        bondDrawingLimit_(bondDrawingLimit) {
    createBonds(atomsVector3D);
}

void Bonds3D::createBonds(AtomsVector3D *atomsVector3D) {
    for (long i = 0; i < atomsVector3D->numberOfEntities(); ++i) {
        for (long j = i+1; j < atomsVector3D->numberOfEntities(); ++j) {

            auto atomDistance = Metrics::distance(
                    atomsVector3D->operator[](i).position(),
                    atomsVector3D->operator[](j).position());
            auto addedGuiRadii = (GuiHelper::radiusFromType(atomsVector3D->operator[](i).type())
                                + GuiHelper::radiusFromType(atomsVector3D->operator[](j).type()));

            if (atomDistance - 0.5*addedGuiRadii < bondDrawingLimit_)
                new Bond3D(this, *atomsVector3D->particles3D_[i], *atomsVector3D->particles3D_[j]);
        }
    }
}
