//
// Created by heuer on 06.12.16.
//

#include <cmath>
#include "Cylinder.h"
#include "Helper.h"

Cylinder::Cylinder(Qt3DCore::QEntity *root,
                   QColor color,
                   const std::pair<QVector3D, QVector3D>& pair,
                   const float radius,
                   const float alpha)
  : Abstract3dObject(root, QColor(), MidPointVector(pair)),
    radius_(radius),
    start_(pair.first),
    end_(pair.second) {

  difference_ = end_ - start_;
  length_ = difference_.length();

  mesh_ = new Qt3DExtras::QCylinderMesh;
  mesh_->setRadius(radius);
  mesh_->setLength(length_);
  mesh_->setRings(10);
  mesh_->setSlices(20);
  material->setAlpha(alpha);

  rotateToOrientation(difference_);

  material->setAmbient(color);

  entity->addComponent(mesh_);
}

Cylinder::Cylinder(const Cylinder &cylinder)
        : radius_(cylinder.getRadius()),
          length_(cylinder.getLength()),
          start_(cylinder.getStart()),
          end_(cylinder.getEnd()),
          difference_(cylinder.difference_)
{
}

void Cylinder::rotateToOrientation(const QVector3D &orientation) {

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
    auto tmp = start_;
    start_ = end_;
    end_ = tmp;
  }
}
