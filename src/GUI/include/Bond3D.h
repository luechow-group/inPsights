//
// Created by Michael Heuer on 23.12.16.
//

#ifndef AMOLQCGUI_BOND3D_H
#define AMOLQCGUI_BOND3D_H

#include "DividedCylinder.h"
#include "Atom3D.h"

class Bond3D : public DividedCylinder {
public:
  Bond3D(const Atom3D& src, const Atom3D& dest);
  //Bond3D(Atom* src, Atom* dest);

private:
  const Atom3D& src_,dest_;
};

/* TODO delete if not needed
#include "Cylinder.h"
#include "Atom.h"
class Bond3D : public Cylinder {
public:
  //Bond3D(Qt3DCore::QEntity* root, Atom* src, Atom* dest);
  Bond3D(Qt3DCore::QEntity* root,const Atom& src, const Atom& dest);
};*/

#endif //AMOLQCGUI_BOND3D_H
