//
// Created by heuer on 06.12.16.
//

#ifndef AMOLQCPP_CYLINDER_H
#define AMOLQCPP_CYLINDER_H

#include <Qt3DExtras/QCylinderMesh>
#include "Abstract3dObject.h"

class Cylinder : public Abstract3dObject {

public:
  Cylinder(Qt3DCore::QEntity *root, QColor color,
           const std::pair<QVector3D, QVector3D>& pair,
           const float radius,
           const float alpha = 1.0f);

  ~Cylinder() {};

  float getRadius() const { return radius_; };

  void setRadius(const float radius) {
      radius_ = radius;
      mesh_->setRadius(radius);
  };

  float length() const { return difference().length(); };
  QVector3D start() const{ return start_; };
  QVector3D end() const{ return end_; };
  QVector3D difference() const{ return end_ - start_; };

private:
  void rotateToOrientation(const QVector3D &orientation);

  float radius_;
  QVector3D start_, end_;
  Qt3DExtras::QCylinderMesh *mesh_;
};

#endif //AMOLQCPP_CYLINDER_H
