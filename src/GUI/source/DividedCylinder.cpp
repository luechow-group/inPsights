//
// Created by Michael Heuer on 24.12.16.
//

#include "DividedCylinder.h"
#include "Helper.h"

DividedCylinder::DividedCylinder(Qt3DCore::QEntity *root,
                                 const std::pair<QColor,QColor>& colorPair,
                                 const std::pair<QVector3D, QVector3D>& locationPair,
                                 const float& radius)
  : srcCylinder(root,colorPair.first,{locationPair.first, MidPointVector(locationPair)}, radius),
    destCylinder(root,colorPair.second,{MidPointVector(locationPair), locationPair.second}, radius)
{
}
