//
// Created by heuer on 06.12.16.
//

#include <iostream>
#include <Abstract3dObject.h>
#include <Eigen/Core>

Abstract3dObject::Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D& location, float alpha)
  : QEntity(root),
    material(new Qt3DExtras::QPhongAlphaMaterial(this)),
    transform(new Qt3DCore::QTransform),
    picker(new Qt3DRender::QObjectPicker),
    color_(color)
{
  material->setSpecular(Qt::white);
  material->setShininess(0);
  material->setAmbient(color);
  material->setAlpha(alpha);
  transform->setTranslation(location);

  addComponent(transform);
  addComponent(material);
  addComponent(picker);

  picker->setHoverEnabled(true);
}

