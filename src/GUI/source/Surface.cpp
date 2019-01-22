//
// Created by Michael Heuer on 2019-01-07.
//

#include <Surface.h>
#include <SurfaceData.h>

Surface::Surface(
        Qt3DCore::QEntity *root,
        const SurfaceData & surfaceData,
        QColor color, float alpha)
        :
        Abstract3dObject(root, std::move(color), {0,0,0}),
        mesh_(new SurfaceMesh(surfaceData)) {

    addComponent(mesh_);
}
