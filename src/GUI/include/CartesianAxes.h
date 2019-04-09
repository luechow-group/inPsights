//
// Created by heuer on 09.04.19.
//

#ifndef INPSIGHTS_CARTESIANAXES_H
#define INPSIGHTS_CARTESIANAXES_H

#include "Abstract3dObject.h"
#include "Arrow.h"

class CartesianAxes : public Abstract3dObject{
public:
    CartesianAxes(Qt3DCore::QEntity *root, QVector3D origin = {0,0,0},
            float length = 1.0f,
            float baseRadius = 0.025f,
            float alpha = 0.25f);

private:
    Arrow *x_, *y_, *z_;
};

#endif //INPSIGHTS_CARTESIANAXES_H
