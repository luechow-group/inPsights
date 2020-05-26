/* Copyright (C) 2019 Michael Heuer.
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


Motif::Motif(std::set<Eigen::Index> electronIndices, MotifType type)
        : type_(type), electronIndices_(std::move(electronIndices)) {};

Motif::Motif(std::set<Eigen::Index> electronIndices,
             std::set<Eigen::Index> atomIndices,
             MotifType type)
        : type_(type), electronIndices_(std::move(electronIndices)), atomIndices_(std::move(atomIndices)) {};

bool Motif::containsElectronQ(Eigen::Index i) const {
    return std::find(electronIndices_.begin(), electronIndices_.end(), i) != electronIndices_.end();
}

bool Motif::containsAtomQ(Eigen::Index i) const {
    return std::find(atomIndices_.begin(), atomIndices_.end(), i) != atomIndices_.end();
}

bool Motif::operator<(const Motif &rhs) const {
    return electronIndices_ < rhs.electronIndices_;
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

const std::set<Eigen::Index> &Motif::electronIndices() const {
    return electronIndices_;
}

void Motif::setElectronIndices(const std::set<Eigen::Index> &electronIndices) {
    electronIndices_ = electronIndices;
}

const std::set<Eigen::Index> &Motif::atomIndices() const {
    return atomIndices_;
}

void Motif::setAtomIndices(const std::set<Eigen::Index> &atomIndices) {
    Motif::atomIndices_ = atomIndices;
}

namespace YAML {
    Node convert<Motif>::encode(const Motif &rhs) {
        Node node;
        node["Type"] = toString(rhs.type());
        node["ElectronIndices"] =
                std::list<Eigen::Index>(std::begin(rhs.electronIndices()),std::end(rhs.electronIndices()));
        node["AtomIndices"] =
                std::list<Eigen::Index>(std::begin(rhs.atomIndices()),std::end(rhs.atomIndices()));
        return node;
    }

    bool convert<Motif>::decode(const Node &node, Motif &rhs) {
        if (!node.IsMap())
            return false;

        auto electronIndices = node["ElectronIndices"].as<std::list<Eigen::Index>>();
        auto atomIndices = node["AtomIndices"].as<std::list<Eigen::Index>>();

        rhs.setElectronIndices(
                std::set<Eigen::Index>(std::begin(electronIndices),std::end(electronIndices)));
        rhs.setAtomIndices(
                std::set<Eigen::Index>(std::begin(atomIndices),std::end(atomIndices)));
        rhs.setType(fromString(node["Type"].as<std::string>()));
        return true;
    }

    Emitter &operator<<(Emitter &out, const Motif &rhs) {
        out << YAML::Flow << BeginMap
            << Key << "ElectronIndices" << Value << BeginSeq;
        for(auto i : rhs.electronIndices())
            out <<  i;
        out << EndSeq << Key << "AtomIndices" << Value << BeginSeq;
        for(auto i : rhs.atomIndices())
            out <<  i;
        out << EndSeq
        << Key << "Type" << Value << toString(rhs.type())
        << EndMap;
        return out;
    };
}
