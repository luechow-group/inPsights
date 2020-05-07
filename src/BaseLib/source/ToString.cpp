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

#include <iomanip>
#include <cmath>
#include "ToString.h"

std::string ToString::longToString(long a,
                                   unsigned leadingSpaces) {
    std::string string = std::to_string(a);
    string = std::string(leadingSpaces + 1 - string.length(),' ') + string;
    return string;
}

std::string ToString::vector3dToString(const Eigen::Vector3d &vector,
                                       unsigned decimalPlaces, unsigned leadingSpaces) {
    Eigen::VectorXd vectorXd(3);
    vectorXd << vector;
    return vectorXdToString(vectorXd,decimalPlaces,leadingSpaces);
}

std::string ToString::doubleToString(double a, unsigned decimalPlaces, unsigned leadingSpaces) {
    std::stringstream sstream;
    if (a >= 0){
        sstream << " ";
    }
    for (unsigned i = 1; i <= leadingSpaces; i++){
        if (fabs(a) < pow(10.0,i)){
            sstream << " ";
        }
    }
    sstream << std::fixed << std::setprecision(decimalPlaces) << a;
    return sstream.str();
}

std::string ToString::vectorXdToString(const Eigen::VectorXd &vector, unsigned decimalPlaces, unsigned leadingSpaces) {
    std::stringstream sstream;
    for (Eigen::Index i = 0; i < vector.size(); i++){
        sstream << " " << doubleToString(vector(i), decimalPlaces, leadingSpaces);
    }
    return sstream.str();
}

std::string ToString::matrixXdToString(const Eigen::MatrixXd &matrix, unsigned decimalPlaces, unsigned leadingSpaces) {
    std::stringstream ss;
    for (Eigen::Index i = 0; i < matrix.rows(); i++) {
        ss << ToString::vectorXdToString(matrix.row(i), decimalPlaces, leadingSpaces) << std::endl;
    }
    return ss.str();
}

std::string ToString::vectorXiToString(const Eigen::VectorXi &vector) {
    std::stringstream ss;
    ss << vector.transpose();
    return ss.str();
}

std::string ToString::stdvectorIntToString(const std::vector<int>& vector) {
    std::stringstream ss;
    for(auto i : vector)
        ss << i << " ";
    return ss.str();
}

std::string ToString::stdvectorUIntToString(const std::vector<unsigned >& vector) {
    std::stringstream ss;
    for(auto i : vector)
        ss << i << " ";
    return ss.str();
}

std::string ToString::stdvectorLongIntToString(const std::vector<long int>& vector) {
    std::stringstream ss;
    for(auto i : vector)
        ss << i << " ";
    return ss.str();
}

std::string ToString::stdvectorLongUIntToString(const std::vector<long unsigned>& vector) {
    std::stringstream ss;
    for(auto i : vector)
        ss << i << " ";
    return ss.str();
}