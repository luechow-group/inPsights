// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CARTESIANAXES_H
#define INPSIGHTS_CARTESIANAXES_H

#include "Abstract3dObject.h"
#include "Arrow.h"

class CartesianAxes : public Abstract3dObject{
public:
    CartesianAxes(Qt3DCore::QEntity *root, QVector3D origin = {0,0,0},
            float length = 5.0f,
            float baseRadius = 0.01f,
            float alpha = 0.75f);

private:
    Arrow *x_, *y_, *z_;
};

#endif //INPSIGHTS_CARTESIANAXES_H
