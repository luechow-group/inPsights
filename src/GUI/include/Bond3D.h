// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_BOND3D_H
#define INPSIGHTS_BOND3D_H

#include "DividedCylinder.h"
#include "Particle3D.h"
#include <Particle.h>

class Bond3D : public DividedCylinder {
public:
  Bond3D(Qt3DCore::QEntity* root, const Atom3D& src, const Atom3D& dest); //TODO why use const ref?

};

#endif //INPSIGHTS_BOND3D_H
