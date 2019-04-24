/* Copyright (C) 2018-2019 Michael Heuer.
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

#ifndef INPSIGHTS_MODIFIEDSPHERICALBESSER1STKIND_H
#define INPSIGHTS_MODIFIEDSPHERICALBESSER1STKIND_H

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

#endif //INPSIGHTS_MODIFIEDSPHERICALBESSER1STKIND_H
