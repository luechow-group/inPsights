//
// Created by heuer on 06.12.16.
//

#include "Sphere.h"

Sphere::Sphere(Qt3DCore::QEntity *root, QColor color, const QVector3D location, const float radius)
  : Abstract3dObject(root, color, location),
    radius_(radius) {

  mesh_ = new Qt3DExtras::QSphereMesh;
  mesh_->setRadius(radius);
  mesh_->setRings(10);
  mesh_->setSlices(20);

  entity->addComponent(mesh_);
}
