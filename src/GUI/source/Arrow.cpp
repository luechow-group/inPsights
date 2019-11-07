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

#include "Arrow.h"
#include <cassert>

Arrow::Arrow(Qt3DCore::QEntity *root,
             QColor color,
             const std::pair<QVector3D, QVector3D> &pair,
             float baseRadius,
             float relativeTipLength,
             float relativeTipBottomRadius,
             float alpha)
        :
        Cylinder(root, color, {pair.first,calculateReducedCylinderEnd(pair, relativeTipLength)},
                 baseRadius, alpha),
        tip_(new Cone(root, color, {calculateReducedCylinderEnd(pair, relativeTipLength), pair.second}, baseRadius * relativeTipBottomRadius, 0.0f, alpha)) {
    assert(baseRadius > 0);
    assert(relativeTipLength > 0 && "The relative tip length must be a positive value.");
    assert(relativeTipLength < 0.5f && "The relative tip length must be shorter than the half length of the vector.");
    assert(relativeTipBottomRadius > 1 && "The relative tip bottom radius must be a number greater than one.");
}

QVector3D Arrow::calculateReducedCylinderEnd(const std::pair<QVector3D, QVector3D> &pair, float relativeTipLength) {
    assert(relativeTipLength > 0 && "The relative tip length must be a positive value.");
    assert(relativeTipLength < 0.5f && "The relative tip length must be shorter than the half length of the vector.");

    auto difference = pair.second - pair.first;
    return pair.first + difference * (1.0f-relativeTipLength);
}
