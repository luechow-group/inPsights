//
// Created by Leonard Reuter on 09.03.18.
//

#ifndef AMOLQCGUI_TOSTRING_H
#define AMOLQCGUI_TOSTRING_H

#include <string>
#include <Eigen/Core>

namespace ToString {
    std::string unsignedLongToString(unsigned long a,
                                     unsigned leadingSpaces = 1);

    std::string doubleToString(double a,
                               unsigned decimalPlaces = 5,
                               unsigned leadingSpaces = 1);

    std::string vector3dToString(const Eigen::Vector3d &vector,
                                 unsigned decimalPlaces = 5,
                                 unsigned leadingSpaces = 1);

    std::string vectorXdToString(const Eigen::VectorXd &vector,
                                 unsigned decimalPlaces = 5,
                                 unsigned leadingSpaces = 1);
}

#endif //AMOLQCGUI_TOSTRING_H
