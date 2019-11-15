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

#include <cmath>
#include "Cylinder.h"
#include "GuiHelper.h"

Cylinder::Cylinder(Qt3DCore::QEntity *root,
                   QColor color,
                   const std::pair<QVector3D, QVector3D>& pair,
                   float radius,
                   float alpha)
  : Abstract3dObject(root, std::move(color), GuiHelper::midPointVector(pair), alpha),
    radius_(radius),
    start_(pair.first),
    end_(pair.second),
    mesh_(new Qt3DExtras::QCylinderMesh(this)) {

    mesh_->setRadius(radius);
    mesh_->setLength(length());
    mesh_->setRings(8);
    mesh_->setSlices(16);

    rotateToOrientation(difference());

    addComponent(mesh_);
}

float Cylinder::length() const {
    return difference().length();
}

QVector3D Cylinder::start() const {
    return start_;
}

QVector3D Cylinder::end() const {
    return end_;
}
QVector3D Cylinder::difference() const{
    return end_ - start_;
};

float Cylinder::getRadius() const {
    return radius_;
};

void Cylinder::setRadius(const float radius) {
    radius_ = radius;
    mesh_->setRadius(radius);
};

void Cylinder::addToXml(std::ostream &os, unsigned int sortKey) const  {
    auto cylinderColor = color();
    auto center = transform->translation();

    QVector3D axis;
    float angle;
    transform->rotation().getAxisAndAngle(&axis,&angle);

    os << "<transform translation='"
       << center[0] << ","
       << center[1] << ","
       << center[2]
       <<"' rotation='"
       << axis[0] << ","
       << axis[1] << ","
       << axis[2] << ","
       << angle*ConversionFactors::deg2rad
       << "'>\n";
    os << "<shape><appearance sortKey='"
       << sortKey
       << "'><material diffuseColor='"
       << cylinderColor.redF() << " "
       << cylinderColor.greenF() << " "
       << cylinderColor.blueF()
       << "' transparency='" << 0.75*(1.0f-material->alpha()) << "'></material></appearance>\n";

    os << "<cylinder height='"
       << length()
       << "' radius='"
       << getRadius()
       << "' top='false' bottom='false'></cylinder></shape></transform>\n\n";
}

void Cylinder::rotateToOrientation(const QVector3D &orientation) {

  auto origVec = QVector3D(0, 1, 0);

  auto perpendicular = QVector3D::crossProduct(origVec, orientation);
  auto l = perpendicular.length();
  float epsilon = 0.00001;

  if (l > epsilon) {
    QVector3D axis = perpendicular.normalized();
    float angle = float(atan2f(l, QVector3D::dotProduct(origVec, orientation)) / (M_PI * 2) * 360.0f);
    transform->setRotation(QQuaternion::fromAxisAndAngle(axis, angle));
  } else if (QVector3D::dotProduct(origVec, orientation) > 0.0f) {
    // Nearly positively aligned; skip rotation, or compute
    // axis and angle using other means
  } else {
    // Nearly negatively aligned; axis is any vector perpendicular
    // to either vector, and angle is 180 degrees
    auto tmp = start_;
    start_ = end_;
    end_ = tmp;
  }
}
