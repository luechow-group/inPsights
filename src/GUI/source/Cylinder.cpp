//
// Created by heuer on 06.12.16.
//

#include <cmath>

#include "Cylinder.h"
#include "Helper.h"

Cylinder::Cylinder(Qt3DCore::QEntity *root,
                   QColor color,
                   const std::pair<QVector3D, QVector3D> pair,
                   const float radius)
  : Abstract3dObject(root, QColor(), MidPointVector(pair)),
    radius(radius),
    start(pair.first),
    end(pair.second) {

  difference = end - start;
  length = difference.length();

  mesh = new Qt3DExtras::QCylinderMesh;
  mesh->setRadius(radius);
  mesh->setLength(length);
  mesh->setRings(100);
  mesh->setSlices(10);

  rotateToOrientation(difference);

  material->setAmbient(color);

  entity->addComponent(mesh);
}

void Cylinder::rotateToOrientation(const QVector3D orientation) {

  auto origVec = QVector3D(0, 1, 0);

  auto perpendicular = QVector3D::crossProduct(origVec, orientation);
  auto l = perpendicular.length();
  float epsilon = 0.00001;

  if (l > epsilon) {
    QVector3D axis = perpendicular.normalized();
    float angle = float(atan2f(l, QVector3D::dotProduct(origVec, orientation)) / (M_PI * 2) * 360.0f);
    transform->setRotation(QQuaternion::fromAxisAndAngle(axis, angle));
  } else if (QVector3D::dotProduct(origVec, orientation) > 0.0f) {
    // Nearly positively aligned; skip rotation, or compute
    // axis and angle using other means
  } else {
    // Nearly negatively aligned; axis is any vector perpendicular
    // to either vector, and angle is 180 degrees
    auto tmp = start;
    start = end;
    end = tmp;
  }
}
