/* Copyright (C) 2016-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
