//
// Created by heuer on 06.12.16.
//

#ifndef AMOLQC_ABSTRACT3DOBJECT_H
#define AMOLQC_ABSTRACT3DOBJECT_H

#include <QVector3D>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QSphereMesh>
#include <Qt3DRender/QObjectPicker>

class Abstract3dObject : public Qt3DCore::QEntity {
  Q_OBJECT
public:
  Abstract3dObject(){};
  Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D location);
  //~Abstract3dObject(){};

  Qt3DCore::QEntity* entity;
  Qt3DExtras::QPhongMaterial* material;
  Qt3DCore::QTransform* transform;
  Qt3DRender::QObjectPicker *picker;

  //void setColor(const QColor& color){ this->color = color;};
  //void setLocation(const QVector3D& location){ this->location = location;};

  QColor getColor() const { return color; };
  QVector3D getLocation() const { return location; };

public slots:
  void onPressed(bool pressed);

protected:
  QColor color;
  QVector3D location;
};

#endif //AMOLQC_ABSTRACT3DOBJECT_H
