// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_LINE3D_H
#define INPSIGHTS_LINE3D_H

#include "Abstract3dObject.h"
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QBuffer>
#include <Qt3DRender/QGeometry>
#include <Qt3DRender/QGeometryRenderer>
#include "GuiHelper.h"

class Line3D : public Abstract3dObject {

public:
    Line3D(Qt3DCore::QEntity *root, QColor color,
           const std::pair<QVector3D, QVector3D> &pair,
           float alpha = 1.0f);
    
    float length() const;

    QVector3D start() const;

    QVector3D end() const;

    QVector3D difference() const;

private:
    QVector3D start_, end_;

    Qt3DRender::QGeometry *createGeometry(Qt3DCore::QEntity *entity, const std::pair<QVector3D, QVector3D> &pair) const;
};

#endif //INPSIGHTS_LINE3D_H
