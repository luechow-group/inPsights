//
// Created by Michael Heuer on 08.05.18.
//

#include "Particle.h"

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