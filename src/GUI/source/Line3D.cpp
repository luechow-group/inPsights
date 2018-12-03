//
// Created by heuer on 03.12.18.
//

#include <Line3D.h>

Line3D::Line3D(Qt3DCore::QEntity *root, QColor color,
               const std::pair<QVector3D, QVector3D> &pair,
               float alpha)
        :
        Abstract3dObject(root, std::move(color), GuiHelper::midPointVector(pair), alpha),
        start_(pair.first),
        end_(pair.second) {

    // mesh
    auto *line = new Qt3DRender::QGeometryRenderer(this);
    line->setGeometry(createGeometry(this, pair));
    line->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);

    addComponent(line);
}

float Line3D::length() const { return difference().length(); };

QVector3D Line3D::start() const { return start_; };

QVector3D Line3D::end() const { return end_; };

QVector3D Line3D::difference() const { return end_ - start_; };

Qt3DRender::QGeometry *
Line3D::createGeometry(Qt3DCore::QEntity *entity, const std::pair<QVector3D, QVector3D> &pair) const {
    auto *geometry = new Qt3DRender::QGeometry(entity);

    // position vertices (start and end)
    QByteArray bufferBytes;
    bufferBytes.resize(3 * 2 * sizeof(float)); // start.x, start.y, start.end + end.x, end.y, end.z
    auto *positions = reinterpret_cast<float *>(bufferBytes.data());

    *positions++ = (pair.first - transform->translation()).x();
    *positions++ = (pair.first - transform->translation()).y();
    *positions++ = (pair.first - transform->translation()).z();
    *positions++ = (pair.second - transform->translation()).x();
    *positions++ = (pair.second - transform->translation()).y();
    *positions++ = (pair.second - transform->translation()).z();

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
    auto *indices = reinterpret_cast<unsigned int *>(indexBytes.data());
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
