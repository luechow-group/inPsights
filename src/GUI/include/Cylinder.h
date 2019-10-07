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

  float getRadius() const { return radius_; };

  void setRadius(const float radius) {
      radius_ = radius;
      mesh_->setRadius(radius);
  };

  float length() const { return difference().length(); };
  QVector3D start() const{ return start_; };
  QVector3D end() const{ return end_; };
  QVector3D difference() const{ return end_ - start_; };


  friend std::ostream& operator<< (std::ostream& os, const Cylinder& obj) {
      auto color = obj.color();
      auto center = obj.transform->translation();

      QVector3D axis;
      float angle;
      obj.transform->rotation().getAxisAndAngle(&axis,&angle);

      os << "<transform translation='"
         << center[0] << ","
         << center[1] << ","
         << center[2]
         <<"' rotation='"
         << axis[0] << ","
         << axis[1] << ","
         << axis[2] << ","
         << angle*ConversionFactors::deg2rad
         << "'>\n";
      os << "<shape><appearance><material diffuseColor='"
         << color.red() << " "
         << color.green() << " "
         << color.blue()
         << "' transparency='" << obj.material->alpha() << "'></material></appearance>\n";

      os << "<cylinder top='false' bottom='false' height='"
         << obj.length()
         << "' radius='"
         << obj.getRadius()
         << "'></cylinder>\n";
      os <<"</shape></transform>\n\n";
      return os;
  }

private:
  void rotateToOrientation(const QVector3D &orientation);

  float radius_;
  QVector3D start_, end_;
  Qt3DExtras::QCylinderMesh *mesh_;
};

#endif //INPSIGHTS_CYLINDER_H
