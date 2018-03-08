//
// Created by Michael Heuer on 07.03.18.
//

#ifndef AMOLQCPP_POSITIONFORMAT_H
#define AMOLQCPP_POSITIONFORMAT_H

#include <Eigen/Core>

namespace PositionFormat{
    static const unsigned significantDigits = 6;
    static const std::string separator = " ";
    static const Eigen::IOFormat positionFormat = Eigen::IOFormat(significantDigits, 0, separator, "\n", "", "", "", "");
};

#endif //AMOLQCPP_POSITIONFORMAT_H
