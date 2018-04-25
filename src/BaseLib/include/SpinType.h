//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_SPINTYPE_H
#define AMOLQCPP_SPINTYPE_H

#include <string>
#include <cassert>

namespace Spins {
    enum class SpinType {
        none=0, alpha=-1, beta=-2
    };

    SpinType first();

    SpinType last();

    SpinType spinTypeFromInt(int type);

    int spinTypeToInt(SpinType spinType);

    std::string toString(const SpinType& s);

    double magneticQuantumNumber(SpinType spinType);

    double quantumNumber();
}

std::ostream& operator<<(std::ostream& os, const Spins::SpinType& s);

#endif //AMOLQCPP_SPINTYPE_H
