//
// Created by heuer on 09.12.16.
//

#include <Atom3D.h>
#include <GuiHelper.h>

Atom3D::Atom3D(Qt3DCore::QEntity *root, const QVector3D& location, const Element& elementType)
  : Sphere(root,
           GuiHelper::QColorFromType<Element>(elementType),
           location,
           GuiHelper::radiusFromType(elementType)),
  elementType_(elementType) {
  material->setAlpha(0.25f);
  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}


Atom3D::Atom3D(const Atom3D &atom)
  : Sphere(atom.parentEntity(),
           atom.getColor(),
           atom.getLocation(),
           atom.getRadius()),
  elementType_(atom.getElementType()){
  material->setAlpha(0.25f);
  //connect(picker, &Qt3DRender::QObjectPicker::pressedChanged, this, &Atom::onPressed);
}

/*void Atom::onPressed(bool pressed) {
  if (pressed) std::cout << "pressed" << std::endl;
  else std::cout << "not pressed" << std::endl;
}*/
