//
// Created by heuer on 06.12.16.
//

#include <iostream>
#include <Abstract3dObject.h>
#include <Eigen/Core>

Abstract3dObject::Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D& location)
  : QEntity(root),
    color_(color)
{
  entity = new Qt3DCore::QEntity(root);
  material = new Qt3DExtras::QPhongAlphaMaterial(root);
  transform = new Qt3DCore::QTransform;
  picker = new Qt3DRender::QObjectPicker;

  material->setSpecular(Qt::white);
  material->setShininess(0);
  material->setAmbient(color);
  material->setAlpha(1.0f);
  transform->setTranslation(location);

  entity->addComponent(transform);
  entity->addComponent(material);
  entity->addComponent(picker);


  picker->setHoverEnabled(true);
}

