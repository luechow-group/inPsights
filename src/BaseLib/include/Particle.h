//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_PARTICLE_H
#define AMOLQCPP_PARTICLE_H

#include <Eigen/Core>
#include <iostream>
#include "ToString.h"
#include "ElementType.h"
#include "ElementInfo.h"
#include "SpinType.h"

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

using Electron = Particle<Spins::SpinType>;
using Atom = Particle<Element>;

template<>
std::string Electron::toString() const;
template<>
std::string Atom::toString() const;

template<>
int Electron::charge() const;
template<>
int Atom::charge() const;


#endif //AMOLQCPP_PARTICLE_H
