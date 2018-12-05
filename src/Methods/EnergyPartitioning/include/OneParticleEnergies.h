//
// Created by Michael Heuer on 22.11.18.
//

#ifndef INPSIGHTS_ONEPARTICLEENERGIES_H
#define INPSIGHTS_ONEPARTICLEENERGIES_H


#include <iostream>
#include <math.h>
#include <yaml-cpp/yaml.h>

namespace OneParticleEnergies {
    void oneAtomEnergies(const YAML::Node &cluster, const YAML::Node &Vnn);

    void oneElectronEnergies(const YAML::Node &cluster);
}

#endif //INPSIGHTS_ONEPARTICLEENERGIES_H
