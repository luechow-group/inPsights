//
// Created by Michael Heuer on 05.02.18.
//

#include "ElementType.h"
#include "ElementInfo.h"

Elements::ElementType Elements::first() {
    return Elements::ElementType::H;
};

Elements::ElementType Elements::last(){
    return Elements::ElementType::Og;
};

Elements::ElementType Elements::elementTypeFromInt(int type){
    return static_cast<Elements::ElementType>(type);
};

int Elements::elementTypeToInt(Elements::ElementType elementType){
    return int(elementType);
};

std::ostream& operator<< (std::ostream& os, const Elements::ElementType & e){
    os <<  Elements::ElementInfo::symbol(e);
    return os;
};
