//
// Created by heuer on 10.05.17.
//

#include "Electron3D.h"

Electron3D::Electron3D(Qt3DCore::QEntity *root, const QVector3D& location, const Spin::SpinType& spinType)
        : Sphere(root,
                 Spin::QColorFromSpinType(spinType),
                 location,
                 float(Elements::ElementInfo::vdwRadius(Elements::ElementType::H)/10.0f/4.0f)),
          spinType_(spinType) {

  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}


Electron3D::Electron3D(const Electron3D &electron)
        : Sphere(electron.parentEntity(),
                 electron.getColor(),
                 electron.getLocation(),
                 electron.getRadius()),
          spinType_(electron.getSpinType()){

  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}