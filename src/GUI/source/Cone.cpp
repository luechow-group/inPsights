// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include "Cone.h"
#include "GuiHelper.h"

Cone::Cone(Qt3DCore::QEntity *root,
           QColor color,
           const std::pair<QVector3D, QVector3D>& pair,
           float bottomRadius,
           float topRadius,
           float alpha)
        : Abstract3dObject(root, QColor(), GuiHelper::midPointVector(pair)),
          topRadius_(topRadius),
          bottomRadius_(bottomRadius),
          start_(pair.first),
          end_(pair.second),
          mesh_(new Qt3DExtras::QConeMesh(this)) {

  mesh_->setTopRadius(topRadius);
  mesh_->hasTopEndcapChanged(true);
  mesh_->setBottomRadius(bottomRadius);
  mesh_->hasBottomEndcapChanged(true);
  mesh_->setLength(length());
  mesh_->setRings(100);
  mesh_->setSlices(10);

  material->setAlpha(alpha);

  rotateToOrientation(difference());

  material->setAmbient(color);

  addComponent(mesh_);
}

void Cone::rotateToOrientation(const QVector3D &orientation) {

  auto origVec = QVector3D(0, 1, 0);

  auto perpendicular = QVector3D::crossProduct(origVec, orientation);
  auto l = perpendicular.length();
  double epsilon = 1.0e-14f;

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
