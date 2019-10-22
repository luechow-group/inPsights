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

#include "Cylinder.h"
#include "Cone.h"

#ifndef INPSIGHTS_POLYLINE_H
#define INPSIGHTS_POLYLINE_H

class Polyline : public Abstract3dObject{

public:
    Polyline(Qt3DCore::QEntity *root, QColor color, const std::vector<QVector3D> points, const float radius,
             bool arrowTipQ = false);
private:
    std::vector<QVector3D> points_;
    std::vector<Cylinder*> cylinders_;
    float radius_;
    float totalArcLength_;
    bool arrowTipQ_;
    Cone* arrowTip_;
};

#endif //INPSIGHTS_POLYLINE_H
