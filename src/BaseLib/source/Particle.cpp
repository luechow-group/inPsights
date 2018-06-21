//
// Created by Michael Heuer on 08.05.18.
//

#include "Particle.h"
#include <yaml-cpp/yaml.h>

template<> std::string Electron::toString() const {
    std::string typeString = "";
    typeString += "e" + Spins::toString(Spins::spinFromInt(type_));
    typeString += ToString::vector3dToString(position_);

    return typeString;
}

template<> std::string Atom::toString() const {
    std::string typeString = "";
    typeString += Elements::ElementInfo::symbol(Elements::elementFromInt(type_));
    if(typeString.length() == 1)
        typeString += " ";
    typeString += ToString::vector3dToString(position_);

    return typeString;
}

template<>
int Electron::charge() const{
    return -1;
}
template<>
int Atom::charge() const{
    return Elements::ElementInfo::Z(Elements::elementFromInt(type_));
}

namespace YAML {

    Node convert<Eigen::Vector3d>::encode(const Eigen::Vector3d &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        return node;
    }

    bool convert<Eigen::Vector3d>::decode(const Node &node, Eigen::Vector3d &rhs) {
        if (!node.IsSequence() || node.size() != 3) {
            return false;
        }
        rhs[0] = node[0].as<double>();
        rhs[1] = node[1].as<double>();
        rhs[2] = node[2].as<double>();
        return true;
    }

    Emitter &operator<<(Emitter &out, const Eigen::Vector3d &v) {
        out << Flow << BeginSeq << v[0] << v[1] << v[2] << EndSeq;
        return out;
    }

    //}namespace YAML {

    //    struct convert<Particle<Type>> {

    Node convert<TypedParticle>::encode(const TypedParticle& rhs) {
        Node node;
        node.push_back(rhs.type());
        node.push_back(rhs.position());
        return node;
    }
    Node convert<Atom>::encode(const Atom& rhs) {
        Node node;
        node.push_back(Elements::ElementInfo::symbol(rhs.type()));
        node.push_back(rhs.position());
        return node;
    }
    Node convert<Electron>::encode(const Electron& rhs) {
        Node node;
        node.push_back(Spins::toString(rhs.type()));
        node.push_back(rhs.position());
        return node;
    }

    bool convert<TypedParticle>::decode(const Node& node, TypedParticle& rhs) {
        if(!node.IsSequence() || node.size() != 2) {
            return false;
        }
        rhs = {node[0].as<int>(),node[0].as<Eigen::Vector3d>()};
        return true;
    }
    bool convert<Atom>::decode(const Node& node, Atom& rhs) {
        if(!node.IsSequence() || node.size() != 2) {
            return false;
        }
        rhs = {node[0].as<Element>(),node[0].as<Eigen::Vector3d>()};
        return true;
    }
    bool convert<Electron>::decode(const Node& node, Electron& rhs) {
        if(!node.IsSequence() || node.size() != 2) {
            return false;
        }
        rhs = {node[0].as<Spin>(),node[0].as<Eigen::Vector3d>()};
        return true;
    }

    Emitter& operator<< (Emitter& out, const TypedParticle & p) {
        out << Flow << BeginSeq << p.type() << p.position() << EndSeq;
        return out;
    }
    Emitter& operator<< (Emitter& out, const Atom& p) {
        out << Flow << BeginSeq << p.type() << p.position() << EndSeq;
        return out;
    }
    Emitter& operator<< (Emitter& out, const Electron & p) {
        out << Flow << BeginSeq << p.type() << p.position() << EndSeq;
        return out;
    }
}