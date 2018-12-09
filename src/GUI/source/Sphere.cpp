//
// Created by heuer on 06.12.16.
//

#include <utility>
#include <Sphere.h>
#include "Sphere.h"

Sphere::Sphere(Qt3DCore::QEntity *root, QColor color, const QVector3D location, const float radius)
        :
        Abstract3dObject(root, std::move(color), location),
        highlightedQ_(false),
        selectedQ_(false),
        radius_(radius),
        mesh_(new Qt3DExtras::QSphereMesh(this)) {

    mesh_->setRadius(radius);
    mesh_->setRings(8);
    mesh_->setSlices(16);

    addComponent(mesh_);

    //QObject::connect(picker, &Qt3DRender::QObjectPicker::containsMouseChanged, this, &Sphere::onHighlighted);
}

void Sphere::onHighlighted(bool highlightQ) {
    highlightedQ_ = highlightQ;
    update();
}

void Sphere::onSelected(bool selectedQ) {
    selectedQ_ = selectedQ;
    update();
}

void Sphere::update() {
    if(highlightedQ_) {
        mesh_->setRadius(radius_ * 1.25f);
        material->setAmbient(QColor(0, 255, 255));
    } else if (selectedQ_) {
        mesh_->setRadius(radius_);
        material->setAmbient(QColor(255, 0, 255));
    } else {
        mesh_->setRadius(radius_);
        material->setAmbient(color());
    }
}
