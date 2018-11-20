//
// Created by Michael Heuer on 24.12.16.
//

#ifndef INPSIGHTS_DIVIDEDCYLINDER_H
#define INPSIGHTS_DIVIDEDCYLINDER_H

#include "Cylinder.h"

class DividedCylinder : Abstract3dObject{
public:
  DividedCylinder(Qt3DCore::QEntity *root,
                  const std::pair<QColor,QColor>& colorPair,
                  const std::pair<QVector3D, QVector3D>& locationPair,
                  float radius,
                  float alpha = 1.0f);

  ~DividedCylinder(){}


private:
  Cylinder srcCylinder_,destCylinder_;
};

#endif //INPSIGHTS_DIVIDEDCYLINDER_H
