// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <ParticleIndices.h>

ParticleIndices::ParticleIndices(const std::set<Eigen::Index> &indices)
: particleIndices_(indices) {}

const std::set<Eigen::Index> &ParticleIndices::indices() const {
    return particleIndices_;
}

std::set<Eigen::Index>& ParticleIndices::indices() {
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

bool ParticleIndices::operator==(const ParticleIndices& rhs) {
    return this->particleIndices_ == rhs.particleIndices_;
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