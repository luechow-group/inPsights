//
// Created by heuer on 06.12.16.
//

#include "Sphere.h"

Sphere::Sphere(Qt3DCore::QEntity *root, QColor color, const QVector3D location, const float radius)
  : Abstract3dObject(root, color, location),
    radius(radius) {

  mesh = new Qt3DExtras::QSphereMesh;
  mesh->setRadius(radius);
  mesh->setRings(100);
  mesh->setSlices(100);

  entity->addComponent(mesh);
}
