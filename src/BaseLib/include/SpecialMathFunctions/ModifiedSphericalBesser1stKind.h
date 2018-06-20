//
// Created by Michael Heuer on 20.06.18.
//

#ifndef AMOLQCPP_MODIFIEDSPHERICALBESSER1STKIND_H
#define AMOLQCPP_MODIFIEDSPHERICALBESSER1STKIND_H

#include <vector>

struct ModifiedSphericalBessel1stKind
{
    ModifiedSphericalBessel1stKind(int degree);
    void evaluate(double r, bool differentiate);

    static std::vector<double> eval(int degree, double r);
    static constexpr double RADZERO = 1e-10;
    static constexpr double SPHZERO = 1e-4;

    int _degree;
    std::vector<double> _in;
    std::vector<double> _din;
};

#endif //AMOLQCPP_MODIFIEDSPHERICALBESSER1STKIND_H
