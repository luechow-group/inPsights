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

#include <iostream>
#include <Abstract3dObject.h>
#include <Eigen/Core>

Abstract3dObject::Abstract3dObject(Qt3DCore::QEntity *root, QColor color, const QVector3D& location, float alpha)
  : QEntity(root),
    material(new Qt3DExtras::QPhongAlphaMaterial(this)),
    transform(new Qt3DCore::QTransform),
    picker(new Qt3DRender::QObjectPicker),
    color_(std::move(color))
{
  material->setSpecular(Qt::white);
  material->setShininess(4.0f);
  material->setAmbient(color);
  material->setAlpha(alpha);
  transform->setTranslation(location);

  addComponent(transform);
  addComponent(material);
  addComponent(picker);

  picker->setHoverEnabled(true);
}

QColor Abstract3dObject::color() const {
    return color_;
}
