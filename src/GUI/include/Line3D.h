//
// Created by Michael Heuer on 19.11.18.
//

#ifndef INPSIGHTS_LINE3D_H
#define INPSIGHTS_LINE3D_H

#include "Abstract3dObject.h"
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QGeometryRenderer>
#include "GuiHelper.h"

class Line3D : public Abstract3dObject {

public:
    Line3D(Qt3DCore::QEntity *root, QColor color,
           const std::pair<QVector3D, QVector3D> &pair,
           float alpha = 1.0f)
            :
            Abstract3dObject(root, std::move(color), GuiHelper::midPointVector(pair), alpha),
            start_(pair.first),
            end_(pair.second) {

        // mesh
        auto *line = new Qt3DRender::QGeometryRenderer(this->parentNode());
        line->setGeometry(getGeometry(pair));
        line->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);

        // entity
        addComponent(line);
    }

    Qt3DRender::QGeometry *getGeometry( const std::pair<QVector3D, QVector3D> &pair) const {
        auto *geometry = new Qt3DRender::QGeometry(this->parentNode());

        // position vertices (start and end)
        QByteArray bufferBytes;
        bufferBytes.resize(3 * 2 * sizeof(float)); // start.x, start.y, start.end + end.x, end.y, end.z
        float *positions = reinterpret_cast<float *>(bufferBytes.data());

        *positions++ = pair.first.x();
        *positions++ = pair.first.y();
        *positions++ = pair.first.z();
        *positions++ = pair.second.x();
        *positions++ = pair.second.y();
        *positions++ = pair.second.z();

        auto *buf = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, geometry);
        buf->setData(bufferBytes);

        auto *positionAttribute = new Qt3DRender::QAttribute(geometry);
        positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());
        positionAttribute->setVertexBaseType(Qt3DRender::QAttribute::Float);
        positionAttribute->setVertexSize(3);
        positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
        positionAttribute->setBuffer(buf);
        positionAttribute->setByteStride(3 * sizeof(float));
        positionAttribute->setCount(2);
        geometry->addAttribute(positionAttribute); // We add the vertices in the geometry

        // connectivity between vertices
        QByteArray indexBytes;
        indexBytes.resize(2 * sizeof(unsigned int)); // start to end
        unsigned int *indices = reinterpret_cast<unsigned int *>(indexBytes.data());
        *indices++ = 0;
        *indices++ = 1;

        auto *indexBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, geometry);
        indexBuffer->setData(indexBytes);

        auto *indexAttribute = new Qt3DRender::QAttribute(geometry);
        indexAttribute->setVertexBaseType(Qt3DRender::QAttribute::UnsignedInt);
        indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
        indexAttribute->setBuffer(indexBuffer);
        indexAttribute->setCount(2);
        geometry->addAttribute(indexAttribute); // We add the indices linking the points in the geometry

        return geometry;
    };


    float length() const { return difference().length(); };

    QVector3D start() const { return start_; };

    QVector3D end() const { return end_; };

    QVector3D difference() const { return end_ - start_; };

private:
    QVector3D start_, end_;
};

#endif //INPSIGHTS_LINE3D_H
