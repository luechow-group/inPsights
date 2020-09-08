// Copyright (C) 2018 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

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