//
// Created by heuer on 06.12.16.
//

#ifndef INPSIGHTS_ABSTRACT3DOBJECT_H
#define INPSIGHTS_ABSTRACT3DOBJECT_H

#include <QVector3D>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DExtras/QPhongAlphaMaterial>
#include <Qt3DExtras/QSphereMesh>
#include <Qt3DRender/QObjectPicker>

class Abstract3dObject : public Qt3DCore::QEntity {
  Q_OBJECT
public:
  Abstract3dObject(){};
  Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D& location);
  //~Abstract3dObject(){};

  Qt3DCore::QEntity* entity;
  Qt3DExtras::QPhongAlphaMaterial* material;
  Qt3DCore::QTransform* transform;
  //Qt3DRender::QObjectPicker *picker;

  //void setColor(const QColor& color){ this->color = color;};
  //void setLocation(const QVector3D& location){ this->location = location;};

  void setAlpha(float alpha);

  QColor getColor() const { return color_; };
  QVector3D getLocation() const { return location_; };

//public slots:
  //void onPressed(bool pressed);

protected:
  QColor color_;
  float alpha_;
  QVector3D location_;
};

#endif //INPSIGHTS_ABSTRACT3DOBJECT_H
