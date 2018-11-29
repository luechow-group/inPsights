//
// Created by heuer on 10.05.17.
//

#include <Electron3D.h>
#include <GuiHelper.h>

Electron3D::Electron3D(Qt3DCore::QEntity *root, const QVector3D& location, const Spin& spinType)
        : Sphere(root,
                 GuiHelper::QColorFromType<Spin>(spinType),
                 location,
                 float(Elements::ElementInfo::vdwRadius(Element::H)/10.0f/4.0f)),
          spinType_(spinType) {
  setAlpha(0.5f);
  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}


Electron3D::Electron3D(const Electron3D &electron)
        : Sphere(electron.parentEntity(),
                 electron.getColor(),
                 electron.getLocation(),
                 electron.getRadius()),
          spinType_(electron.getSpinType())
{
  setAlpha(0.5f);
  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}