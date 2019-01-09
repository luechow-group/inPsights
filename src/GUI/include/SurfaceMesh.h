//
// Created by Michael Heuer on 2019-01-03.
//

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
