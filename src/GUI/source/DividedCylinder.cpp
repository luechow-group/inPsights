// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "DividedCylinder.h"
#include "GuiHelper.h"

DividedCylinder::DividedCylinder(Qt3DCore::QEntity *root,
                                 const std::pair<QColor,QColor>& colorPair,
                                 const std::pair<QVector3D, QVector3D>& locationPair,
                                 float radius,
                                 float alpha)
  : QEntity(root),
    srcCylinder_(new Cylinder(
            this,
            colorPair.first,
            {locationPair.first, GuiHelper::midPointVector(locationPair)},
            radius,
            alpha)),
    destCylinder_(new Cylinder(
            this,
            colorPair.second,
            {GuiHelper::midPointVector(locationPair),
             locationPair.second},
             radius,
             alpha))
{}
