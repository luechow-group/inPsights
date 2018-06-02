//
// Created by Michael Heuer on 01.02.18.
//
#include "SpinType.h"

Spin Spins::first() {
    return Spin::alpha;
};

Spin Spins::last(){
    return Spin::beta;
};

Spin Spins::spinFromInt(int type){
    return static_cast<Spin>(type);//-storageShift);
};

int Spins::spinToInt(Spin spinType){
    return int(spinType);//+Spin::storageShift;
};

std::string Spins::toString(const Spin& s){
    switch(s) {
        case Spin::alpha: return "a";
        case Spin::beta: return "b";
        case Spin::none: return "-";
    }
}

std::ostream& operator<<(std::ostream& os, const Spin& s){
    os << Spins::toString(s);
    return os;
}

double Spins::magneticQuantumNumber(Spin spinType) {
    assert(spinType != Spin::none && "The Spin::SpinType cannot be 'none'.");

    switch(spinType) {
        case Spin::alpha:
            return 1/2.;
        case Spin::beta:
            return -1/2.;
        default:
            return 0;
    }
};

double Spins::quantumNumber(){
    return 1/2.0;
};