/* Copyright (C) 2019 Michael Heuer.
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

#ifndef INPSIGHTS_ARROW_H
#define INPSIGHTS_ARROW_H

#include "Cylinder.h"
#include "Cone.h"
#include <NaturalConstants.h>

class Arrow : public Cylinder {

public:
    Arrow(Qt3DCore::QEntity *root, QColor color,
             const std::pair<QVector3D, QVector3D>& pair,
             float baseRadius,
             float relativeTipLength = 0.25f,
             float relativeTipBottomRadius = 1.25f,
             float alpha = 1.0f);
private:
    Cone* tip_;

    QVector3D calculateReducedCylinderEnd(const std::pair<QVector3D, QVector3D> &pair, float relativeTipLength);
};

#endif //INPSIGHTS_ARROW_H
