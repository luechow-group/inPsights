//
// Created by heuer on 09.04.19.
//

#include "CartesianAxes.h"

CartesianAxes::CartesianAxes(Qt3DCore::QEntity *root, QVector3D origin, float length, float baseRadius, float alpha)
        :
        Abstract3dObject(root, QColor(), origin),
        x_(new Arrow(this, Qt::red,  {origin - QVector3D(length/2, 0, 0), origin + QVector3D(length/2, 0, 0)}, baseRadius, 0.1f, 2.0f, alpha)),
        y_(new Arrow(this, Qt::green,{origin - QVector3D(0, length/2, 0), origin + QVector3D(0, length/2, 0)}, baseRadius, 0.1f, 2.0f, alpha)),
        z_(new Arrow(this, Qt::blue, {origin - QVector3D(0, 0, length/2), origin + QVector3D(0, 0, length/2)}, baseRadius, 0.1f, 2.0f, alpha))
{}
