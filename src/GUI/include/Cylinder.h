//
// Created by heuer on 06.12.16.
//

#ifndef AMOLQCGUI_CYLINDER_H
#define AMOLQCGUI_CYLINDER_H

#include <Qt3DExtras/QCylinderMesh>
#include "Abstract3dObject.h"

#include "iostream"

class Cylinder : public Abstract3dObject {

public:
  Cylinder(Qt3DCore::QEntity *root, QColor color,
           const std::pair<QVector3D, QVector3D> pair,
           const float radius);

  ~Cylinder() {};

  float getRadius() const { return radius; };
  float getLength() const { return length; };

private:
  void rotateToOrientation(const QVector3D orientation);

  float radius, length;
  QVector3D start, end;
  QVector3D difference;
  Qt3DExtras::QCylinderMesh *mesh;

};

#endif //AMOLQCGUI_CYLINDER_H
