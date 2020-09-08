// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
