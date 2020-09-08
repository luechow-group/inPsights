// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
