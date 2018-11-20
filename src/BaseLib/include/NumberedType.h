//
// Created by Michael Heuer on 04.06.18.
//

#ifndef INPSIGHTS_NUMBEREDTYPE_H
#define INPSIGHTS_NUMBEREDTYPE_H

#include "ElementType.h"
#include "SpinType.h"
#include <iostream>

template <typename Type>
class NumberedType{
public:
    NumberedType() = default;
    NumberedType(Type type, unsigned number)
            : type_(type),number_(number) {}

    bool operator==(NumberedType<Type> other) const {
        return (type_ == other.type_) && (number_ == other.number_);
    }

    bool operator<(const NumberedType<Type>& other) const {
        if (type_ < other.type_)
            return true;
        else if (type_ == other.type_)
            return (number_ < other.number_);
        else
            return false;
    }

    NumberedType<int> toIntType(){
        return {int(type_),number_};
    }

    friend std::ostream& operator<<(std::ostream& os, const NumberedType& numberedType){
        os << numberedType.number_ << numberedType.type_ << std::endl;
        return os;
    }

    Type type_;
    unsigned number_;
};

using NumberedElement = NumberedType<Element >;
using NumberedSpin = NumberedType<Spin>;

#endif //INPSIGHTS_NUMBEREDTYPE_H
