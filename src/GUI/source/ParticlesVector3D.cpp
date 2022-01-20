// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <ParticlesVector3D.h>
#include <Bonds3D.h>

template <>
void AtomsVector3D::drawConnections(const double &limit) {
    new Bonds3D(this, limit);
}

//void ElectronsVector3D::drawCorrelations() {}
