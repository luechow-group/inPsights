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

#include <Surface.h>
#include <SurfaceData.h>

Surface::Surface(
        Qt3DCore::QEntity *root,
        const SurfaceData & surfaceData,
        QColor color, float alpha)
        :
        Abstract3dObject(root, std::move(color), {0,0,0}, alpha),
        mesh_(new SurfaceMesh(surfaceData)) {

    addComponent(mesh_);
}
