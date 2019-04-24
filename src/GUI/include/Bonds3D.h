/* Copyright (C) 2018-2019 Michael Heuer.
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

#ifndef INPSIGHTS_BONDS3D_H
#define INPSIGHTS_BONDS3D_H

#include <IConnection.h>
#include <ParticlesVector3D.h>
#include <NaturalConstants.h>

class Bonds3D : public IConnection {
public:
    explicit Bonds3D(AtomsVector3D *atomsVector3D, double bondDrawingLimit = 1.40 * 1e-10 / AU::length);

private:
    double bondDrawingLimit_;
    void createBonds(AtomsVector3D *atomsVector3D);
};

#endif //INPSIGHTS_BONDS3D_H
