//
// Created by heuer on 09.12.16.
//

#ifndef AMOLQCGUI_ATOM3D_H
#define AMOLQCGUI_ATOM3D_H

#include "Sphere.h"
#include "ElementInfo.h"

class Atom3D : public Sphere {
  //Q_OBJECT
public:
  Atom3D(const Atom3D& atom);
  Atom3D(Qt3DCore::QEntity *root,
       const QVector3D& location,
       const Elements::ElementType& elementType);

  Elements::ElementType getElementType() const { return elementType_; };

  Qt3DRender::QObjectPicker *picker;

//public slots:
  //void onPressed(bool pressed);

private:
  const Elements::ElementType elementType_;
};

#endif //AMOLQCGUI_ATOM3D_H
