/* Copyright (C) 2017-2019 Michael Heuer.
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

#ifndef INPSIGHTS_SPINTYPE_H
#define INPSIGHTS_SPINTYPE_H

#include <string>
#include <cassert>

namespace Spins {
    enum class SpinType {
        alpha=-2, beta=-1,none=0 // use negative numbers to avoid clash with elments
    };

    SpinType first();

    SpinType last();

    SpinType spinFromInt(int type);

    int spinToInt(SpinType spinType);

    std::string toString(const SpinType& s);

    Spins::SpinType fromString(const std::string& s);

    double magneticQuantumNumber(SpinType spinType);

    double quantumNumber();
}

using Spin = Spins::SpinType;

std::ostream& operator<<(std::ostream& os, const Spins::SpinType& s);

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<Spin> {
        static Node encode(const Spin &rhs);
        static bool decode(const Node &node, Spin &rhs);
    };
    Emitter &operator<<(Emitter &out, const Spin &e);
}

#endif //INPSIGHTS_SPINTYPE_H
