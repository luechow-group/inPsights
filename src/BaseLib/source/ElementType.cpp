//
// Created by Michael Heuer on 05.02.18.
//

#include "ElementType.h"
#include "ElementInfo.h"

Element Elements::first() {
    return Element::H;
};

Element Elements::last(){
    return Element::Og;
};

Element Elements::elementFromInt(int type){
    return static_cast<Element>(type);
};

int Elements::elementToInt(Element element){
    return int(element);
};

std::ostream& operator<< (std::ostream& os, const Element & e){
    os <<  Elements::ElementInfo::symbol(e);
    return os;
};
