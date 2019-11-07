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

#include <ParticlesVector3D.h>
#include <Bonds3D.h>
#include <SameSpinConnections3D.h>

template <>
void AtomsVector3D::drawConnections() {
    new Bonds3D(this);
}

template <>
void ElectronsVector3D::drawConnections() {
    new SameSpinConnections3D(this);
}

//void ElectronsVector3D::drawCorrelations() {}