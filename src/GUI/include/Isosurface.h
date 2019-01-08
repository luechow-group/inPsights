//
// Created by Michael Heuer on 2019-01-07.
//

#ifndef INPSIGHTS_ISOSURFACE_H
#define INPSIGHTS_ISOSURFACE_H

#include <Abstract3dObject.h>
#include <IsosurfaceMesh.h>

class Isosurface : public Abstract3dObject{
public:
    Isosurface(Qt3DCore::QEntity *root,
               const std::vector<Vertex>& vertices,
               const std::vector<Triangle>& triangles,
            QColor color, float alpha = 1.0 );

private:
    IsosurfaceMesh * mesh_;
};

#endif //INPSIGHTS_ISOSURFACE_H
