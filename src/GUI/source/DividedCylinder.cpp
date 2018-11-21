//
// Created by Michael Heuer on 24.12.16.
//

#include "DividedCylinder.h"
#include "GuiHelper.h"

DividedCylinder::DividedCylinder(Qt3DCore::QEntity *root,
                                 const std::pair<QColor,QColor>& colorPair,
                                 const std::pair<QVector3D, QVector3D>& locationPair,
                                 float radius,
                                 float alpha)
  : srcCylinder_(root,colorPair.first,{locationPair.first, GuiHelper::midPointVector(locationPair)}, radius, alpha),
    destCylinder_(root,colorPair.second,{GuiHelper::midPointVector(locationPair), locationPair.second}, radius, alpha)
{
}
