//
// Created by Michael Heuer on 23.12.16.
//

#ifndef INPSIGHTS_BOND3D_H
#define INPSIGHTS_BOND3D_H

#include "DividedCylinder.h"
#include "Particle3D.h"
#include <Particle.h>

class Bond3D : public DividedCylinder {
public:
  Bond3D(const Atom3D& src, const Atom3D& dest);

private:
};

#endif //INPSIGHTS_BOND3D_H
