// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ENUMERATEDTYPE_H
#define INPSIGHTS_ENUMERATEDTYPE_H

#include "ElementType.h"
#include "SpinType.h"
#include <iostream>

template <typename Type>
class EnumeratedType{
public:
    EnumeratedType() = default;
    EnumeratedType(Type type, unsigned number)
            : type_(type),number_(number) {}

    bool operator==(EnumeratedType<Type> other) const {
        return (type_ == other.type_) && (number_ == other.number_);
    }

    bool operator<(const EnumeratedType<Type>& other) const {
        if (type_ < other.type_)
            return true;
        else if (type_ == other.type_)
            return (number_ < other.number_);
        else
            return false;
    }

    EnumeratedType<int> toIntType(){
        return {int(type_),number_};
    }

    friend std::ostream& operator<<(std::ostream& os, const EnumeratedType& enumeratedType){
        os << enumeratedType.number_ << enumeratedType.type_ << std::endl;
        return os;
    }

    Type type_;
    unsigned number_;
};

using EnumeratedElement = EnumeratedType<Element >;
using EnumeratedSpin = EnumeratedType<Spin>;

#endif //INPSIGHTS_ENUMERATEDTYPE_H
