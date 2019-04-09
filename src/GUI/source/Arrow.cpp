//
// Created by heuer on 09.04.19.
//

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
