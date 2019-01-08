//
// Created by Michael Heuer on 2019-01-03.
//

#ifndef INPSIGHTS_ISOSURFACEMESH_H
#define INPSIGHTS_ISOSURFACEMESH_H

#include <Qt3DRender/QGeometryRenderer>
#include <vector>
#include <Vertex.h>
#include <Triangle.h>


class IsosurfaceMesh : public Qt3DRender::QGeometryRenderer {
public:
    explicit IsosurfaceMesh(const std::vector<Vertex>& vertices,
                            const std::vector<Triangle>& triangles,
                            Qt3DCore::QNode *parent = nullptr);

    const unsigned short nCoordinates = 3; // cartesian coordinates
    const unsigned short nIndicesPerTriangle= 3;
};



#endif //INPSIGHTS_ISOSURFACEMESH_H
