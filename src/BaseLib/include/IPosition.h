// Copyright (C) 2019 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later


#ifndef INPSIGHTS_IPOSITIONS_H
#define INPSIGHTS_IPOSITIONS_H

#include <ParticlesVector.h>

namespace YAML{
    Eigen::Vector3d decodePosition(const YAML::Node &node, const AtomsVector &nuclei);
}

#endif //INPSIGHTS_IPOSITIONS_H
