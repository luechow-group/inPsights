// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
