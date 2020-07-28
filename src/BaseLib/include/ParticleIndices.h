/* Copyright (C) 2020 Michael Heuer.
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
