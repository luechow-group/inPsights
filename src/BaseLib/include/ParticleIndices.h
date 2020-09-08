// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <set>
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>
#include "ElementType.h"
#include "SpinType.h"

#ifndef INPSIGHTS_PARTICLEGROUPINDICES_H
#define INPSIGHTS_PARTICLEGROUPINDICES_H

class ParticleIndices{
public:
    ParticleIndices() = default;
    ParticleIndices(const std::set<Eigen::Index>& indices);

    const std::set<Eigen::Index> &indices() const;
    std::set<Eigen::Index> &indices();

    void setIndices(const std::set<Eigen::Index>& indices);

    void merge(const ParticleIndices& other);

    bool containsIndexQ(Eigen::Index i) const;

    bool operator==(const ParticleIndices& rhs);

private:
    std::set<Eigen::Index> particleIndices_;
};

namespace YAML {
    class Node; class Emitter;
    template<> struct convert<ParticleIndices> {
        static Node encode(const ParticleIndices &rhs);
        static bool decode(const Node &node, ParticleIndices &rhs);
    };
    Emitter &operator<<(Emitter &out, const ParticleIndices &rhs) ;
}

#endif //INPSIGHTS_PARTICLEGROUPINDICES_H
