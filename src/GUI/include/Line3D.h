/* Copyright (C) 2018-2019 Michael Heuer.
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
