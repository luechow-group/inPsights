// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SurfaceMesh.h"
#include <SurfaceData.h>

#include <Qt3DCore/QEntity>
#include <Qt3DRender/QRenderAspect>
#include <Qt3DRender/QGeometryRenderer>
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QBuffer>

SurfaceMesh::SurfaceMesh(
        const SurfaceData& surfaceData,
        Qt3DCore::QNode *parent) {

    // Custom Mesh
    auto *customGeometry = new Qt3DRender::QGeometry(this);

    auto *vertexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, customGeometry);
    auto *indexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, customGeometry);

    // Vertices
    auto vertexDataPackageSize = nCoordinates*2; // position + normal vector
    QByteArray vertexBufferData;
    vertexBufferData.resize(surfaceData.vertices.size() * vertexDataPackageSize * sizeof(float));
    auto *rawVertexArray = reinterpret_cast<float *>(vertexBufferData.data());

    int idx = 0;
    for (const auto & v : surfaceData.vertices) {
        rawVertexArray[idx++] = float(v.position[0]);
        rawVertexArray[idx++] = float(v.position[1]);
        rawVertexArray[idx++] = float(v.position[2]);

        rawVertexArray[idx++] = float(v.normal[0]);
        rawVertexArray[idx++] = float(v.normal[1]);
        rawVertexArray[idx++] = float(v.normal[2]);
    }
    vertexDataBuffer->setData(vertexBufferData);

    auto *positionAttribute = new Qt3DRender::QAttribute();
    positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    positionAttribute->setBuffer(vertexDataBuffer);
    positionAttribute->setDataType(Qt3DRender::QAttribute::Float);
    positionAttribute->setDataSize(nCoordinates);
    positionAttribute->setByteOffset(0);
    positionAttribute->setByteStride(vertexDataPackageSize  * sizeof(float));
    positionAttribute->setCount(surfaceData.vertices.size());
    positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());

    auto *normalAttribute = new Qt3DRender::QAttribute();
    normalAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    normalAttribute->setBuffer(vertexDataBuffer);
    normalAttribute->setDataType(Qt3DRender::QAttribute::Float);
    normalAttribute->setDataSize(nCoordinates);
    normalAttribute->setByteOffset(nCoordinates * sizeof(float));
    normalAttribute->setByteStride(vertexDataPackageSize * sizeof(float));
    normalAttribute->setCount(surfaceData.vertices.size());
    normalAttribute->setName(Qt3DRender::QAttribute::defaultNormalAttributeName());

    QByteArray indexBufferData;
    indexBufferData.resize(surfaceData.triangles.size() * nIndicesPerTriangle * sizeof(unsigned));
    auto *rawIndexArray = reinterpret_cast<unsigned*>(indexBufferData.data());
    idx = 0;
    for (const auto& t : surfaceData.triangles) {

        rawIndexArray[idx++] = t.indices[0];
        rawIndexArray[idx++] = t.indices[1];
        rawIndexArray[idx++] = t.indices[2];
        //std::cout << t.vertices[0].index <<", "<< t.vertices[1].index <<", "<< t.vertices[2].index << std::endl;
    }
    indexDataBuffer->setData(indexBufferData);

    auto *indexAttribute = new Qt3DRender::QAttribute();
    indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
    indexAttribute->setBuffer(indexDataBuffer);
    indexAttribute->setDataType(Qt3DRender::QAttribute::UnsignedInt);
    indexAttribute->setDataSize(1);
    indexAttribute->setByteOffset(0);
    indexAttribute->setByteStride(0);
    indexAttribute->setCount(surfaceData.triangles.size()*nIndicesPerTriangle);

    customGeometry->addAttribute(positionAttribute);
    customGeometry->addAttribute(normalAttribute);
    customGeometry->addAttribute(indexAttribute);

    setInstanceCount(1);
    setFirstVertex(0);
    setIndexOffset(0);
    setFirstInstance(0);
    setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
    setGeometry(customGeometry);
    setVertexCount(surfaceData.triangles.size()*nIndicesPerTriangle);
}
