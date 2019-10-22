/* Copyright (C) 2016-2019 Michael Heuer.
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

#ifndef INPSIGHTS_ABSTRACT3DOBJECT_H
#define INPSIGHTS_ABSTRACT3DOBJECT_H

#include <QVector3D>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DExtras/QPhongAlphaMaterial>
#include <Qt3DRender/QObjectPicker>
#include <QToolTip>

class Abstract3dObject : public Qt3DCore::QEntity {
public:
    Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D &location, float alpha = 1.0);

    QColor color() const;

    Qt3DExtras::QPhongAlphaMaterial *material;
    Qt3DCore::QTransform *transform;
public:
    Qt3DRender::QObjectPicker *picker;

private:
    QColor color_;
};

#endif //INPSIGHTS_ABSTRACT3DOBJECT_H
