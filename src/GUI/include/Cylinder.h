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
