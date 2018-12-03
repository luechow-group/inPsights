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
  : QEntity(root),
    srcCylinder_(new Cylinder(
            static_cast<Qt3DCore::QEntity*>(this),
            colorPair.first,
            {locationPair.first, GuiHelper::midPointVector(locationPair)},
            radius,
            alpha)),
    destCylinder_(new Cylinder(
            static_cast<Qt3DCore::QEntity*>(this),
            colorPair.second,
            {GuiHelper::midPointVector(locationPair),
             locationPair.second},
             radius,
             alpha))
{}
