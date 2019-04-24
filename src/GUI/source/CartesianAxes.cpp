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

#include "CartesianAxes.h"

CartesianAxes::CartesianAxes(Qt3DCore::QEntity *root, QVector3D origin, float length, float baseRadius, float alpha)
        :
        Abstract3dObject(root, QColor(), origin),
        x_(new Arrow(this, Qt::red,  {origin - QVector3D(length/2, 0, 0), origin + QVector3D(length/2, 0, 0)}, baseRadius, 0.1f, 2.0f, alpha)),
        y_(new Arrow(this, Qt::green,{origin - QVector3D(0, length/2, 0), origin + QVector3D(0, length/2, 0)}, baseRadius, 0.1f, 2.0f, alpha)),
        z_(new Arrow(this, Qt::blue, {origin - QVector3D(0, 0, length/2), origin + QVector3D(0, 0, length/2)}, baseRadius, 0.1f, 2.0f, alpha))
{}
