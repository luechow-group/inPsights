//
// Created by Michael Heuer on 04.06.18.
//

#ifndef AMOLQCPP_NUMBEREDTYPE_H
#define AMOLQCPP_NUMBEREDTYPE_H

template <typename Type>
class NumberedType{
public:
    NumberedType(Type type, unsigned number)
            : type_(type),number_(number) {}

    bool operator==(NumberedType<Type> other) const {
        return (type_ == other.type_) && (number_ == other.number_);
    }

    bool operator<(const NumberedType<Type>& other) const {
        return (type_ < other.type_) || (number_ < other.number_);
    }

    NumberedType<int> toIntType(){
        return {int(type_),number_};
    }

    friend std::ostream& operator<<(std::ostream& os, const NumberedType& numberedType){
        os << numberedType.number_ << numberedType.type_ << std::endl;
        return os;
    }

    const Type type_;
    const unsigned number_;
};

using NumberedElement = NumberedType<Element >;
using NumberedSpin = NumberedType<Spin>;

#endif //AMOLQCPP_NUMBEREDTYPE_H
