//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_SPINTYPE_H
#define AMOLQCPP_SPINTYPE_H

#include <string>
#include <cassert>

namespace Spins {
    enum class SpinType {
        alpha=-2, beta=-1,none=0
    };

    SpinType first();

    SpinType last();

    SpinType spinTypeFromInt(int type);

    int spinTypeToInt(SpinType spinType);

    std::string toString(const SpinType& s);

    double magneticQuantumNumber(SpinType spinType);

    double quantumNumber();
}

using Spin = Spins::SpinType;

std::ostream& operator<<(std::ostream& os, const Spins::SpinType& s);

#endif //AMOLQCPP_SPINTYPE_H
