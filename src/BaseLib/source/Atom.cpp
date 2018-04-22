//
// Created by heuer on 24.05.17.
//

#include "Atom.h"

#include "ElementInfo.h"
#include "ToString.h"

using namespace Eigen;
using namespace Elements;

Atom::Atom(const Vector3d& position, const ElementType& elementType)
        : Particle(position,ElementInfo::Z(elementType)) {}

ElementType Atom::elementType() const {
    return static_cast<ElementType>(type_);
};

void Atom::setElementType(const ElementType &elementType) {
    type_ = ElementInfo::Z(elementType);
}

int Atom::charge() const {
    return type_;
}

std::ostream& operator<< (std::ostream& os, const Atom& atom) {
    os << atom.toString();
    return os;
}

std::string Atom::toString() const {
    std::string string;
    if (elementType() != ElementType::none){
        string = ElementInfo::symbol(elementType());
        if(string.length() == 1){
            string += " ";
        }
    }
    else{
        string = "  ";
    }

    return string + ToString::vector3dToString(position_);
}