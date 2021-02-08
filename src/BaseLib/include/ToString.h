// Copyright (C) 2018 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_TOSTRING_H
#define INPSIGHTS_TOSTRING_H

#include <string>
#include <Eigen/Core>

namespace ToString {
    std::string longToString(long a,
                             unsigned leadingSpaces = 1);

    std::string doubleToString(double a,
                               unsigned decimalPlaces = 5,
                               unsigned leadingSpaces = 1,
                               bool spaceForPositiveNumber = true);

    std::string vector3dToString(const Eigen::Vector3d &vector,
                                 unsigned decimalPlaces = 5,
                                 unsigned leadingSpaces = 1);

    std::string vectorXdToString(const Eigen::VectorXd &vector,
                                 unsigned decimalPlaces = 5,
                                 unsigned leadingSpaces = 1);

    std::string vectorXiToString(const Eigen::VectorXi &vector);

    std::string matrixXdToString(const Eigen::MatrixXd &matrix,
                                 unsigned decimalPlaces = 5,
                                 unsigned leadingSpaces = 1);

    std::string stdvectorIntToString(const std::vector<int>& vector);
    std::string stdvectorUIntToString(const std::vector<unsigned>& vector);
    std::string stdvectorLongIntToString(const std::vector<long int>& vector);
    std::string stdvectorLongUIntToString(const std::vector<long unsigned>& vector);
}

#endif //INPSIGHTS_TOSTRING_H
