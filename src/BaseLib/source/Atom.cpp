//
// Created by heuer on 24.05.17.
//

#include "Atom.h"
#include "ElementInfo.h"
#include "ToString.h"

using namespace Eigen;
using namespace Elements;

Atom::Atom(const Particle &particle, const ElementType& elementType)
        : Particle(particle),
          elementType_(elementType) {}

Elements::ElementType Atom::elementType() const {
    return elementType_;
};

void Atom::setElementType(const ElementType &elementType) {
    elementType_ = elementType;
}

int Atom::charge() const {
    return int(ElementInfo::Z(elementType_));
}

std::ostream& operator<< (std::ostream& os, const Atom& atom) {
    os << atom.toString();
    return os;
}

std::string Atom::toString() const {
    std::string string;
    if (elementType_ != ElementType::none){
        string = ElementInfo::symbol(elementType_);
        if(string.length() == 1){
            string += " ";
        }
    }
    else{
        string = "  ";
    }

    return string + ToString::vector3dToString(position_);
}