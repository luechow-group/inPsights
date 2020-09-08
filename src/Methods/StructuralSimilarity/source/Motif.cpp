// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Motif.h"
#include <ParticlesVector.h>
#include <Metrics.h>
#include <set>
#include <string>
#include <utility>
#include <yaml-cpp/yaml.h>

std::string toString(MotifType type) {
    switch(type) {
        case MotifType::Core :
            return "Core";
        case MotifType::Valence :
            return "Valence";
        case MotifType::CoreValence :
            return "CoreValence";
        default:
            return "unassigned";
    }
}

MotifType fromString(const std::string& string) {
    if(string == "Core")
        return MotifType::Core;
    else if (string == "Valence")
        return MotifType::Valence;
    else if (string == "CoreValence")
        return MotifType::CoreValence;
    else
        return MotifType::unassigned;
}

Motif::Motif(ParticleIndices electronIndices, MotifType type)
    :
        MolecularSelection(std::move(electronIndices)),
        type_(type)
{};

Motif::Motif(ParticleIndices electronIndices, ParticleIndices nucleiIndices, MotifType type)
    :
        MolecularSelection(std::move(electronIndices), std::move(nucleiIndices)),
        type_(type)
    {};

bool Motif::operator<(const Motif &rhs) const {
    return electrons_.indices() < rhs.electrons_.indices();
}

bool Motif::operator>(const Motif &rhs) const {
    return rhs < *this;
}

bool Motif::operator<=(const Motif &rhs) const {
    return !(rhs < *this);
}

bool Motif::operator>=(const Motif &rhs) const {
    return !(*this < rhs);
}

MotifType Motif::type() const {
    return type_;
}

void Motif::setType(MotifType type) {
    Motif::type_ = type;
}

namespace YAML {
    Node convert<Motif>::encode(const Motif &rhs) {
        Node node;
        node["Type"] = toString(rhs.type());
        node["ElectronIndices"] =
                std::list<Eigen::Index>(std::begin(rhs.electrons_.indices()),std::end(rhs.electrons_.indices()));
        node["AtomIndices"] =
                std::list<Eigen::Index>(std::begin(rhs.nuclei_.indices()),std::end(rhs.nuclei_.indices()));
        return node;
    }

    bool convert<Motif>::decode(const Node &node, Motif &rhs) {
        if (!node.IsMap())
            return false;

        auto electronIndices = node["ElectronIndices"].as<std::list<Eigen::Index>>();
        auto atomIndices = node["AtomIndices"].as<std::list<Eigen::Index>>();

        rhs.electrons_.setIndices(
                std::set<Eigen::Index>(std::begin(electronIndices),std::end(electronIndices)));
        rhs.nuclei_.setIndices(
                std::set<Eigen::Index>(std::begin(atomIndices),std::end(atomIndices)));
        rhs.setType(fromString(node["Type"].as<std::string>()));
        return true;
    }

    Emitter &operator<<(Emitter &out, const Motif &rhs) {
        out << YAML::Flow << BeginMap
            << Key << "ElectronIndices" << Value << BeginSeq;
        for(auto i : rhs.electrons_.indices())
            out <<  i;
        out << EndSeq << Key << "AtomIndices" << Value << BeginSeq;
        for(auto i : rhs.nuclei_.indices())
            out <<  i;
        out << EndSeq
        << Key << "Type" << Value << toString(rhs.type())
        << EndMap;
        return out;
    };
}
