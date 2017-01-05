//
// Created by Michael Heuer on 23.12.16.
//

#ifndef AMOLQCGUI_BOND_H
#define AMOLQCGUI_BOND_H

#include "DividedCylinder.h"
#include "Atom.h"

class Bond : public DividedCylinder {
public:
  Bond(const Atom& src, const Atom& dest);
  //Bond(Atom* src, Atom* dest);

private:
  const Atom& src,dest;
};

/* TODO delete if not needed
#include "Cylinder.h"
#include "Atom.h"
class Bond : public Cylinder {
public:
  //Bond(Qt3DCore::QEntity* root, Atom* src, Atom* dest);
  Bond(Qt3DCore::QEntity* root,const Atom& src, const Atom& dest);
};*/



#endif //AMOLQCGUI_BOND_H
