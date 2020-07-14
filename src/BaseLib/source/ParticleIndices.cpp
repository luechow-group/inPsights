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

#include <ParticleIndices.h>

ParticleIndices::ParticleIndices(const std::set<Eigen::Index> &indices)
: particleIndices_(indices) {}

const std::set<Eigen::Index> &ParticleIndices::indices() const {
    return particleIndices_;
}

void ParticleIndices::setIndices(const std::set<Eigen::Index> &indices) {
    particleIndices_ = indices;
}

bool ParticleIndices::containsIndexQ(Eigen::Index i) const {
    return std::find(particleIndices_.begin(), particleIndices_.end(), i) != particleIndices_.end();
}

void ParticleIndices::merge(const ParticleIndices &other) {
    particleIndices_.insert(std::begin(other.indices()), std::end(other.indices()));
}


namespace YAML {
    Node convert<ParticleIndices>::encode(const ParticleIndices &rhs) {
        Node node;
        node["Indices"] =
                std::list<Eigen::Index>(std::begin(rhs.indices()),std::end(rhs.indices()));
        return node;
    }

    bool convert<ParticleIndices>::decode(const Node &node, ParticleIndices &rhs) {
        if (!node.IsMap())
            return false;

        auto indices = node["Indices"].as<std::list<Eigen::Index>>();
        rhs.setIndices(std::set<Eigen::Index>(std::begin(indices), std::end(indices)));
        return true;
    }

    Emitter &operator<<(Emitter &out, const ParticleIndices &rhs) {
        out << YAML::Flow << BeginMap << Key << "Indices" << Value << BeginSeq;
        for(auto i : rhs.indices())
            out <<  i;
        out << EndSeq << EndMap;
        return out;
    };
}