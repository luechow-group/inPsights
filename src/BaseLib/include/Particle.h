//
// Created by Michael Heuer on 29.10.17.
//

#ifndef INPSIGHTS_PARTICLE_H
#define INPSIGHTS_PARTICLE_H

#include <Eigen/Core>
#include <iostream>

#include "ToString.h"
#include "ElementType.h"
#include "SpinType.h"
#include "EigenYamlConversion.h"

#include <memory>

template <typename Type> class LinkedParticle{
public:

    LinkedParticle(
            const Eigen::Ref<Eigen::Vector3d> &positionRef,
            int* typeNameRef)
    : position_(positionRef), type_(typeNameRef)
    {}

    Eigen::Ref<Eigen::Vector3d>& positionRef(){
        return position_;
    }

    Eigen::Vector3d position() const {
        return Eigen::Vector3d(position_);
    }

    virtual void setPosition(const Eigen::Vector3d & position){
        position_ = position;
    }

    Type type() const {
        return static_cast<Type>(*type_);
    }

    int* typeRef(){
        return type_;
    }

    void setType(const Type & type){
        *type_ = static_cast<int>(type);
    }

private:
    Eigen::Ref<Eigen::Vector3d> position_;
    int* type_;
};


template<typename Type>
class Particle {
public:
    Particle()
            : Particle(Eigen::Vector3d::Zero(), 0)
    {}

    Particle(Type type, Eigen::Vector3d position)
            : Particle(position, int(type))
    {}

    Particle(Eigen::Vector3d position, int typeId = 0)
            : position_(std::move(position)), type_(typeId)
    {}

    const Eigen::Vector3d& position() const{
        return position_;
    }

    virtual void setPosition(const Eigen::Vector3d& position) {
        position_ = position;
    }

    void translate(const Eigen::Vector3d& shift){
        setPosition(position()+shift);
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

    //Generic particles are not charged
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

#endif //INPSIGHTS_PARTICLE_H
