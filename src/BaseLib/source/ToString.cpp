//
// Created by Leonard Reuter on 09.03.18.
//

#include <iomanip>
#include "ToString.h"

std::string ToString::int2string(const int a,
                                 const unsigned leadingSpaces) {
    std::string string = std::to_string(a);
    string = std::string(leadingSpaces + 1 - string.length(),' ') + string;
    return string;
}

std::string ToString::vector3d2string(const Eigen::Vector3d &vector,
                            const unsigned decimalPlaces, const unsigned leadingSpaces) {
    std::ostringstream sstream;
    for (int i = 0; i < vector.size(); i++){
        sstream << " " << double2string(vector(i),decimalPlaces,leadingSpaces);
    }
    return sstream.str();
}

std::string ToString::double2string(double a, const unsigned decimalPlaces, const unsigned leadingSpaces) {
    std::ostringstream sstream;
    if (a >= 0){
        sstream << " ";
    }
    for (int i = 1; i <= leadingSpaces; i++){
        if (abs(a) < pow(10.0,i)){
            sstream << " ";
        }
    }
    sstream << std::fixed << std::setprecision(decimalPlaces) << a;
    return sstream.str();
}
