//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_SPINTYPE_H
#define AMOLQCGUI_SPINTYPE_H

#include <string>

namespace Spin {
    enum class SpinType { alpha=1, none=0, beta=-1};

    std::string toString(const SpinType& s);

    /*double magneticQuantumNumber(const SpinType& spinType){
        return double(spinType)/2.0;
    };

    double quantumNumber(){
        return 1/2;
    };*/
}

std::ostream& operator<<(std::ostream& os, const Spin::SpinType& s);

#endif //AMOLQCGUI_SPINTYPE_H
