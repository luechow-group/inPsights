//
// Created by heuer on 06.12.16.
//

#ifndef INPSIGHTS_ABSTRACT3DOBJECT_H
#define INPSIGHTS_ABSTRACT3DOBJECT_H

#include <QVector3D>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DExtras/QPhongAlphaMaterial>
#include <Qt3DRender/QObjectPicker>
#include <QToolTip>

class Abstract3dObject : public Qt3DCore::QEntity {
Q_OBJECT
public:
    Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D &location);

    QColor color() const { return color_; };

protected:
    Qt3DCore::QEntity *entity;
    Qt3DExtras::QPhongAlphaMaterial *material;
    Qt3DCore::QTransform *transform;
public:
    Qt3DRender::QObjectPicker *picker;


private:
    QColor color_;
};

#endif //INPSIGHTS_ABSTRACT3DOBJECT_H
