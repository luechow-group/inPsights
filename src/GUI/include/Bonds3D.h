// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
