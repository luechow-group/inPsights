//
// Created by Michael Heuer on 05.02.18.
//

#include "ElementType.h"
#include "ElementInfo.h"

std::ostream& operator<< (std::ostream& os, const Elements::ElementType & e){
    os <<  Elements::ElementInfo::symbol(e);
    return os;
}
