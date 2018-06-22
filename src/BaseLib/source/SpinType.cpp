//
// Created by Michael Heuer on 01.02.18.
//
#include "SpinType.h"
#include <yaml-cpp/yaml.h>

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

Spins::SpinType Spins::fromString(const std::string &s) {
    if(s == "a")
        return Spin::alpha;
    else if(s == "b")
        return Spin::beta;
    else
        return Spin::none;
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


namespace YAML {
    Node convert<Spin>::encode(const Spin &rhs) {
        Node node = YAML::convert<std::string>::encode(Spins::toString(rhs));
        return node;
    }

    bool convert<Spin>::decode(const Node &node, Spin &rhs) {
        if (!node.IsScalar()) {
            return false;
        }
        rhs = Spins::fromString(node.as<std::string>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const Spin &s) {
        out << Spins::toString(s);
        return out;
    }


}
