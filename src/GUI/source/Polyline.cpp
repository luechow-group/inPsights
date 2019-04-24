/* Copyright (C) 2017-2019 Michael Heuer.
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

#include "Polyline.h"
#include "GuiHelper.h"

Polyline::Polyline(Qt3DCore::QEntity *root, QColor color, const std::vector<QVector3D> points, const float radius,
                   bool arrowTipQ)
        : Abstract3dObject(root,color, QVector3D({0,0,0})),
          points_(points),
          cylinders_(0),
          radius_(radius),
          arrowTipQ_(arrowTipQ)
{

  totalArcLength_ = 0;
  for (size_t i = 1; i < points.size(); ++i) {
    totalArcLength_ += (points[i]-points[i-1]).length();
    cylinders_.emplace_back(new Cylinder(this,color,std::make_pair(points[i-1],points[i]),radius));
  }

  if (arrowTipQ) {
    QVector3D arrowTipEnd = points.back();
    QVector3D previousPoint = points[points.size()-2];
    QVector3D arrowTipStart = arrowTipEnd-((arrowTipEnd-previousPoint).normalized()*radius*10.0);

    // Don't plot arrow if the arrow length would be greater than the half total polyline length
    if( (arrowTipEnd-arrowTipStart).length() < totalArcLength_*0.5 )
      arrowTip_ = new Cone(this, color, std::make_pair(arrowTipStart, arrowTipEnd), radius* 4.0f, radius);
  }
}


