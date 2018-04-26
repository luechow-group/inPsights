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

    Particle(const Eigen::Vector3d& position, Type type)
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

    std::string toString() const{
        std::string typeString;

        if(type_ < 0)
            // particle is an electron
            typeString = "e" + Spins::toString(Spins::spinTypeFromInt(type_));
        else if (type_ > 0) {
            // particle is an atom
            typeString = Elements::ElementInfo::symbol(Elements::elementTypeFromInt(type_));
            if(typeString.length() == 1){
                typeString += " ";
            }
        }
        else typeString = "  "; // particle is none

        return typeString + ToString::vector3dToString(position_);
    }

    friend std::ostream& operator<< (std::ostream& os, const Particle<Type>& p) {
        os << p.toString();
        return os;
    }

    int charge() const{
        if(type_ < 0)
            return -1;
        else if (type_ > 0)
            return Elements::ElementInfo::Z(Elements::elementTypeFromInt(type_));
        else
            return 0;
    }
    
protected:
    Eigen::Vector3d position_;
    int type_;
};

using Electron = Particle<Spins::SpinType>;
using Atom = Particle<Elements::ElementType>;

#endif //AMOLQCPP_PARTICLE_H
