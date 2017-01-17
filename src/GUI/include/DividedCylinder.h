//
// Created by Michael Heuer on 24.12.16.
//

#ifndef AMOLQCGUI_DIVIDEDCYLINDER_H
#define AMOLQCGUI_DIVIDEDCYLINDER_H

#include "Cylinder.h"

class DividedCylinder : Abstract3dObject{
public:
  DividedCylinder(Qt3DCore::QEntity *root,
                  const std::pair<QColor,QColor>& colorPair,
                  const std::pair<QVector3D, QVector3D>& locationPair,
                  const float& radius);

  ~DividedCylinder(){}


private:
  Cylinder srcCylinder_,destCylinder_;
};

#endif //AMOLQCGUI_DIVIDEDCYLINDER_H
