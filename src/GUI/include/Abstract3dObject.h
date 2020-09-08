// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
