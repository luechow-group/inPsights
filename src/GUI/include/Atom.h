//
// Created by heuer on 09.12.16.
//

#ifndef AMOLQCGUI_ATOM_H
#define AMOLQCGUI_ATOM_H

#include "Sphere.h"
#include "ElementInfo.h"

class Atom : public Sphere {
  //Q_OBJECT
public:
  Atom(const Atom& atom);
  Atom(Qt3DCore::QEntity *root,
       const QVector3D& location,
       const Elements::ElementType& elementType);

  Elements::ElementType getElementType() const { return elementType; };

  Qt3DRender::QObjectPicker *picker;

//public slots:
  //void onPressed(bool pressed);

private:
  const Elements::ElementType elementType;
};

#endif //AMOLQCGUI_ATOM_H
