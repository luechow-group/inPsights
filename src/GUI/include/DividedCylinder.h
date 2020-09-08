// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_DIVIDEDCYLINDER_H
#define INPSIGHTS_DIVIDEDCYLINDER_H

#include <Qt3DCore/QEntity>
#include "Cylinder.h"

class DividedCylinder : public Qt3DCore::QEntity{
public:
  DividedCylinder(Qt3DCore::QEntity *root,
                  const std::pair<QColor,QColor>& colorPair,
                  const std::pair<QVector3D, QVector3D>& locationPair,
                  float radius,
                  float alpha = 1.0f);

  Cylinder *srcCylinder_, *destCylinder_;
};

#endif //INPSIGHTS_DIVIDEDCYLINDER_H
