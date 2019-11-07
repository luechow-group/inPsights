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
