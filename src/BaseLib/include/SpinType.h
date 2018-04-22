//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_SPINTYPE_H
#define AMOLQCPP_SPINTYPE_H

#include <string>
#include <cassert>

namespace Spin {
    static const int storageShift = -2;

    enum class SpinType { alpha=1, none=0, beta=-1}; // storage shifted: -1,-2,-3

    std::string toString(const SpinType& s);

    static double magneticQuantumNumber(SpinType spinType){
        assert(spinType != SpinType::none && "The SpinType cannot be 'none'.");
        return double(spinType)/2.;
    };

    static double quantumNumber(){
        return 1/2.0;
    };
}

std::ostream& operator<<(std::ostream& os, const Spin::SpinType& s);

#endif //AMOLQCPP_SPINTYPE_H
