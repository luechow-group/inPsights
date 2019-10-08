//
// Created by heuer on 06.12.16.
//

#ifndef INPSIGHTS_CYLINDER_H
#define INPSIGHTS_CYLINDER_H

#include <Qt3DExtras/QCylinderMesh>
#include "Abstract3dObject.h"
#include <iostream>
#include "NaturalConstants.h"

class Cylinder : public Abstract3dObject {

public:
  Cylinder(Qt3DCore::QEntity *root, QColor color,
           const std::pair<QVector3D, QVector3D>& pair,
           float radius,
           float alpha = 1.0f);

  float getRadius() const;
  void setRadius(const float radius);
  float length() const;
  QVector3D start() const;
  QVector3D end() const;
  QVector3D difference() const;

   void addToXml (std::ostream& os, unsigned sortKey = 1) const;

private:
  void rotateToOrientation(const QVector3D &orientation);

  float radius_;
  QVector3D start_, end_;
  Qt3DExtras::QCylinderMesh *mesh_;
};

#endif //INPSIGHTS_CYLINDER_H
