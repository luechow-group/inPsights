#include <utility>
#include <Sphere.h>


//
// Created by heuer on 06.12.16.
//

#include "Sphere.h"

Sphere::Sphere(Qt3DCore::QEntity *root, QColor color, const QVector3D location, const float radius)
        :
        Abstract3dObject(root, std::move(color), location),
        radius_(radius),
        mesh_(new Qt3DExtras::QSphereMesh(this)) {

    mesh_->setRadius(radius);
    mesh_->setRings(8);
    mesh_->setSlices(16);

    addComponent(mesh_);

    QObject::connect(picker, &Qt3DRender::QObjectPicker::containsMouseChanged, this, &Sphere::highlight);
}

void Sphere::highlight(bool highlightQ) {
    if(highlightQ) {
        oldAlpha_ = material->alpha();
        material->setAlpha(1.0f);
    } else {
        material->setAlpha(oldAlpha_);
    }

}
