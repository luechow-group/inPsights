//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_PARTICLE_H
#define AMOLQCPP_PARTICLE_H

#include <Eigen/Core>
#include <iostream>

#include "ToString.h"
#include "ElementType.h"
#include "SpinType.h"
#include "EigenYamlConversion.h"

template<typename Type>
class Particle {
public:
    Particle()
            : position_(Eigen::Vector3d::Zero()),
              type_(0)
    {}

    Particle(Type type, const Eigen::Vector3d &position)
            : position_(position), type_(int(type))
    {}

    Particle(const Eigen::Vector3d& position, int typeId = 0)
            : position_(position), type_(typeId){}

    const Eigen::Vector3d& position() const{
        return position_;
    }

    Eigen::Vector3d& position() {
        return position_;
    }

    Type type() const{
        return static_cast<Type>(type_);
    }

    // Generic implementation
    std::string toString() const{
        std::string typeString = std::to_string(type_);

        if(typeString.size() == 1)
            typeString = "0"+typeString;
        typeString += ToString::vector3dToString(position_);

       return typeString;
    }

    friend std::ostream& operator<< (std::ostream& os, const Particle<Type>& p) {
        os << p.toString();
        return os;
    }

    //Generic
    int charge() const{
        return 0;
    }
    
protected:
    Eigen::Vector3d position_;
    int type_;
};

using TypedParticle = Particle<int>;
using Electron = Particle<Spin>;
using Atom = Particle<Element>;

template<>
std::string Electron::toString() const;
template<>
std::string Atom::toString() const;

template<>
int Electron::charge() const;
template<>
int Atom::charge() const;


namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<TypedParticle> {
        static Node encode(const TypedParticle & rhs);
        static bool decode(const Node& node, TypedParticle & rhs);
    };
    template<> struct convert<Atom> {
        static Node encode(const Atom & rhs);
        static bool decode(const Node& node, Atom& rhs);
    };
    template<> struct convert<Electron> {
        static Node encode(const Electron & rhs);
        static bool decode(const Node& node, Electron& rhs);
    };

    Emitter& operator<< (Emitter& out, const TypedParticle& p);
    Emitter& operator<< (Emitter& out, const Atom& p);
    Emitter& operator<< (Emitter& out, const Electron& p);
}

#endif //AMOLQCPP_PARTICLE_H
