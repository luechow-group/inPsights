//
// Created by Michael Heuer on 02.12.18.
//

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
