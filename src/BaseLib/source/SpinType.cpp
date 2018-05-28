//
// Created by Michael Heuer on 01.02.18.
//
#include "SpinType.h"

Spins::SpinType Spins::first() {
    return Spins::SpinType::alpha;
};

Spins::SpinType Spins::last(){
    return Spins::SpinType::beta;
};

Spins::SpinType Spins::spinTypeFromInt(int type){
    return static_cast<Spins::SpinType>(type);//-storageShift);
};

int Spins::spinTypeToInt(Spins::SpinType spinType){
    return int(spinType);//+Spin::storageShift;
};

std::string Spins::toString(const Spins::SpinType& s){
    switch(s) {
        case Spins::SpinType::alpha: return "a";
        case Spins::SpinType::beta: return "b";
        case Spins::SpinType::none: return "-";
    }
}

std::ostream& operator<<(std::ostream& os, const Spins::SpinType& s){
    os << Spins::toString(s);
    return os;
}

double Spins::magneticQuantumNumber(Spins::SpinType spinType) {
    assert(spinType != Spins::SpinType::none && "The Spin::SpinType cannot be 'none'.");

    switch(spinType) {
        case Spins::SpinType::alpha:
            return 1/2.;
        case Spins::SpinType::beta:
            return -1/2.;
        default:
            return 0;
    }
};

double Spins::quantumNumber(){
    return 1/2.0;
};