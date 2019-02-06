//
// Created by Michael Heuer on 2019-02-06.
//

#include "Motif.h"
#include <ParticlesVector.h>
#include <Metrics.h>
#include <list>
#include <string>
#include <yaml-cpp/yaml.h>

std::string toString(MotifType type) {
    switch(type) {
        case MotifType::Core :
            return "Core";
        case MotifType::Valence :
            return "Valence";
        default:
            return "unassigned";
    }
}

MotifType fromString(const std::string& string) {
    if(string == "Core")
        return MotifType::Core;
    else if (string == "Valence")
        return MotifType::Valence;
    else
        return MotifType::unassigned;
}


Motif::Motif(const std::list<Eigen::Index> &electronIndices, MotifType type)
        : type_(type), electronIndices_(electronIndices) {};

bool Motif::containsQ(Eigen::Index i) const {
    return std::find(electronIndices_.begin(), electronIndices_.end(), i) != electronIndices_.end();
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

void Motif::setType(MotifType type_) {
    Motif::type_ = type_;
}

const std::list<Eigen::Index> &Motif::electronIndices() const {
    return electronIndices_;
}

void Motif::setElectronIndices(const std::list<Eigen::Index> &electronIndices) {
    electronIndices_ = electronIndices;
}

namespace YAML {
    Node convert<Motif>::encode(const Motif &rhs) {
        Node node;
        node["Type"] = toString(rhs.type());
        node["Indices"] = rhs.electronIndices();
        return node;
    }

    bool convert<Motif>::decode(const Node &node, Motif &rhs) {
        if (!node.IsSequence())
            return false;

        rhs.setType(fromString(node["Type"].as<std::string>()));
        rhs.setElectronIndices(node["Indices"].as<std::list<Eigen::Index>>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const Motif &rhs) {
        out << YAML::Flow << BeginMap
            << Key << "Type" << Value << toString(rhs.type())
            << Key << "Indices" << Value << BeginSeq;
        for(auto i : rhs.electronIndices())
            out <<  i;
        auto copy = rhs;
        out << EndSeq << EndMap;
        return out;
    };
}
