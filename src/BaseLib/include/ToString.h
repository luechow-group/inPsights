//
// Created by Leonard Reuter on 09.03.18.
//

#ifndef AMOLQCGUI_TOSTRING_H
#define AMOLQCGUI_TOSTRING_H

#include <string>
#include <Eigen/Core>

namespace ToString {
    std::string int2string(int a,
                           unsigned leadingSpaces = 1);

    std::string double2string(double a,
                              unsigned decimalPlaces = 5,
                              unsigned leadingSpaces = 1);

    std::string vector3d2string(const Eigen::Vector3d &vector,
                                unsigned decimalPlaces = 5,
                                unsigned leadingSpaces = 1);
}

#endif //AMOLQCGUI_TOSTRING_H
