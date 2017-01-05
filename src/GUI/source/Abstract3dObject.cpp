//
// Created by heuer on 06.12.16.
//

#include "Abstract3dObject.h"

#include <iostream>

Abstract3dObject::Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D location)
  : QEntity(root),
    color(color),
    location(location) {

  entity = new Qt3DCore::QEntity(root);
  material = new Qt3DExtras::QPhongMaterial(root);
  transform = new Qt3DCore::QTransform;
  picker = new Qt3DRender::QObjectPicker;

  material->setSpecular(Qt::white);
  material->setShininess(10);
  material->setAmbient(color);
  transform->setTranslation(location);

  entity->addComponent(transform);
  entity->addComponent(material);
  entity->addComponent(picker);

  connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Abstract3dObject::onPressed);
}

void Abstract3dObject::onPressed(bool pressed) {
  if (pressed) std::cout << "pressed" << std::endl;
  else std::cout << "not pressed" << std::endl;
}
