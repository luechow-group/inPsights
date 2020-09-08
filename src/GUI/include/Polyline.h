// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
