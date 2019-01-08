//
// Created by Michael Heuer on 2019-01-07.
//

#include <Isosurface.h>

Isosurface::Isosurface(
        Qt3DCore::QEntity *root,
        const std::vector<Vertex> &vertices,
        const std::vector<Triangle> &triangles,
        QColor color, float alpha)
        :
        Abstract3dObject(root, std::move(color), {0,0,0}),
        mesh_(new IsosurfaceMesh(vertices, triangles)) {

    addComponent(mesh_);
}
