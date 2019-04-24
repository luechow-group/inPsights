/* Copyright (C) 2016-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
