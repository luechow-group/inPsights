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

  float getRadius() const { return radius_; };
  float getLength() const { return length_; };

private:
  void rotateToOrientation(const QVector3D orientation);

  float radius_, length_;
  QVector3D start_, end_;
  QVector3D difference_;
  Qt3DExtras::QCylinderMesh *mesh_;

};

#endif //AMOLQCGUI_CYLINDER_H
