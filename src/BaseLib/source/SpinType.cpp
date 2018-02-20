//
// Created by Michael Heuer on 01.02.18.
//
#include "SpinType.h"

std::string Spin::toString(const Spin::SpinType& s){
    switch(s) {
        case Spin::SpinType::alpha: return "a";
        case Spin::SpinType::beta: return "b";
        case Spin::SpinType::none: return "-";
    }
}

std::ostream& operator<<(std::ostream& os, const Spin::SpinType& s){
    os << Spin::toString(s);
    return os;
}
