//
// Created by leonard on 27.06.19.
//


#ifndef INPSIGHTS_IPOSITIONS_H
#define INPSIGHTS_IPOSITIONS_H

#include <ParticlesVector.h>

namespace YAML{
    Eigen::Vector3d decodePosition(const YAML::Node &node, const AtomsVector &nuclei);
}

#endif //INPSIGHTS_IPOSITIONS_H
