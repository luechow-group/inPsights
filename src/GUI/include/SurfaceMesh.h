/* Copyright (C) 2019 Michael Heuer.
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

#ifndef INPSIGHTS_SURFACEMESH_H
#define INPSIGHTS_SURFACEMESH_H

#include <Qt3DRender/QGeometryRenderer>

class SurfaceData;

class SurfaceMesh : public Qt3DRender::QGeometryRenderer {
public:
    explicit SurfaceMesh(const SurfaceData& surfaceData, Qt3DCore::QNode *parent = nullptr);

    const unsigned short nCoordinates = 3; // cartesian coordinates
    const unsigned short nIndicesPerTriangle= 3;
};



#endif //INPSIGHTS_SURFACEMESH_H
