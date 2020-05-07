/* Copyright (C) 2018 Leonard Reuter.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_TOSTRING_H
#define INPSIGHTS_TOSTRING_H

#include <string>
#include <Eigen/Core>

namespace ToString {
    std::string longToString(long a,
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
