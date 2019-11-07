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

#include <utility>
#include <Sphere.h>

Sphere::Sphere(Qt3DCore::QEntity *root, QColor color, const QVector3D location, const float radius, const float alpha)
        :
        Abstract3dObject(root, std::move(color), location, alpha),
        highlightedQ_(false),
        selectedQ_(false),
        radius_(radius),
        mesh_(new Qt3DExtras::QSphereMesh(this)) {

    mesh_->setRadius(radius);
    mesh_->setRings(16);
    mesh_->setSlices(32);

    addComponent(mesh_);

    //QObject::connect(picker, &Qt3DRender::QObjectPicker::containsMouseChanged, this, &Sphere::onHighlighted);
}

float Sphere::getRadius() const {
    return radius_;
};

void Sphere::setRadius(const float radius) {
    radius_ = radius;
    mesh_->setRadius(radius);
};

void Sphere::onHighlighted(bool highlightQ) {
    highlightedQ_ = highlightQ;
    update();
}

void Sphere::onSelected(bool selectedQ) {
    selectedQ_ = selectedQ;
    update();
}

void Sphere::addToXml(std::ostream &os, unsigned int sortKey) const {
    auto sphereColor = color();
    auto center = transform->translation();

    os << "<transform translation='"
       << center[0] << ","
       << center[1] << ","
       << center[2] << "'>\n";
    os << "<shape><appearance sortKey='"
       << sortKey
       << "'><material diffuseColor='"
       << sphereColor.redF()<< " "
       << sphereColor.greenF() << " "
       << sphereColor.blueF()
       << "' transparency='" << 0.75*(1.0f-material->alpha()) << "'></material></appearance>\n";

    os << "<sphere radius='"
       << getRadius()
       << "'></sphere></shape></transform>\n\n";
}

void Sphere::update() {
    if(highlightedQ_) {
        mesh_->setRadius(radius_ * 1.25f);
        material->setAmbient(QColor(0, 255, 255));
    } else if (selectedQ_) {
        mesh_->setRadius(radius_);
        material->setAmbient(QColor(255, 255, 0));
    } else {
        mesh_->setRadius(radius_);
        material->setAmbient(color());
    }
}
