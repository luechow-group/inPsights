//
// Created by heuer on 10.05.17.
//

#include "Electron.h"

Electron::Electron(Qt3DCore::QEntity *root, const QVector3D& location, const Spin::SpinType& spinType)
        : Sphere(root,
                 Spin::QColorFromSpinType(spinType),
                 location,
                 float(Elements::ElementInfo::vdwRadius(Elements::ElementType::H)/5.)),
          spinType_(spinType) {

  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}


Electron::Electron(const Electron &electron)
        : Sphere(electron.parentEntity(),
                 electron.getColor(),
                 electron.getLocation(),
                 electron.getRadius()),
          spinType_(electron.getSpinType()){

  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}