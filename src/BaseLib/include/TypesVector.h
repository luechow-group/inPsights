//
// Created by Michael Heuer on 22.04.18.
//

#ifndef AMOLQCPP_TYPESVECTOR_H
#define AMOLQCPP_TYPESVECTOR_H

#include <vector>
#include "AbstractVector.h"
#include "SpinType.h"
#include "ElementType.h"
#include "NumberedType.h"
#include "TypesVector.h"
#include "Interval.h"
#include "ReturnAndReset.h"
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>
#include <cassert>

#include "ISliceable.h"

template <typename Type>
class TypesVector : public ISliceable<int> {
public:

    // unsigned long dummy is there to be specialized in TypesVector<Spin::SpinTypes>
    TypesVector(unsigned long size = 0, unsigned long dummy = 0)
    : ISliceable<int>(size)
    {}

    TypesVector(std::vector<Type> types)
    : ISliceable<int>(0)
    {
        for (const auto& type : types){
            assert(int(type) >= int(Spins::first()));
            assert(int(type) <= int(Elements::last()));
            this->append(type);
        }
        resetRef();
    }

    // TEMPALTE IN ISLICABLE?
    TypesVector<Type>& slice(const Interval& interval, const Reset& resetType = Reset::Automatic) {
        ISliceable<int>::slice(interval,resetType); // template specialization necessary?
        return *this;
    }
    TypesVector<Type>& entity(long i, const Reset& resetType = Reset::Automatic) {
        return slice(Interval(i), resetType);
    }

    Type operator[](long i) const {
        return static_cast<Type>(data_[calculateIndex(i)]);
    }

    // TODO MOVE TO ISLICEABLE
    void insert(Type type, long i) {
        assert(i >= 0 && "The index must be positive.");
        assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

        Eigen::VectorXi before = data_.head(i);
        Eigen::VectorXi after = data_.tail(numberOfEntities()-i);

        data_.resize(numberOfEntities()+1);
        data_ << before, int(type), after;

        incrementNumberOfEntities();
        resetRef();
    }
    void prepend(Type type) {
        insert(type,0);
    }
    void append(Type type) {
        insert(type,numberOfEntities());
    }


    bool operator==(const TypesVector<Type> &other) const {
        return ISliceable<int>::operator==(other);
    }
    bool operator!=(const TypesVector<Type> &other) const {
        return ISliceable<int>::operator!=(other);
    }


    unsigned countOccurence(const Type &type) const { // cannot handle slices
        unsigned count = 0;
        for (int i = 0; i < this->numberOfEntities(); ++i) {
            if (this->operator[](i) == type)
                count++;
        }
        return count;
    }

    NumberedType<Type> getNumberedTypeByIndex(long i) const { // cannot handle slices
        auto type = this->operator[](i);
        unsigned count = 0;

        for (long j = 0; j < calculateIndex(i); ++j)
            if(this->operator[](j) == type)
                count++;

        return {type,count};
    };

    std::pair<bool,long> findIndexOfNumberedType(NumberedType<Type> indexedType) const { // cannot handle slices
        assert(indexedType.number_ < numberOfEntities() &&
               "This index is out of bounds.");

        unsigned count = 0;
        for (unsigned j = 0; j < numberOfEntities(); ++j) {
            if(this->operator[](j) == indexedType.type_) count++;
            if(count-1 == indexedType.number_) return {true, j};
        }
        return {false,0};
    }

    //TODO dirty - refactor
    std::vector<std::pair<Type,unsigned>> countTypes() const {

        std::vector<std::pair<Type,unsigned>> typeCountsPair;
        if(numberOfEntities() > 0) {
            if(numberOfEntities() == 1) {
                return {{this->operator[](0),1}};
            }
            auto copy = data_;
            std::sort(copy.data(), copy.data()+numberOfEntities());

            unsigned count = 1;
            for (int i = 1; i < numberOfEntities(); ++i) {
                const auto& lastType = copy[i - 1];
                const auto& currentType = copy[i];

                if (currentType == lastType)
                    count++;
                else {
                    typeCountsPair.push_back({Type(lastType), count}); // cast redundant?
                    count = 1;
                }
                // treat last element
                if(i == numberOfEntities()-1)
                    typeCountsPair.push_back({Type(currentType), count});
            }
        }
        return typeCountsPair;
    }

    friend std::ostream& operator<<(std::ostream& os, const TypesVector& tv){
        for (unsigned long i = 0; i < tv.numberOfEntities(); i++) {
            os << tv[i] << std::endl;
        }
        return os;
    }
};

using IntegerTypesVector = TypesVector<int>;
using SpinTypesVector = TypesVector<Spin>;
using ElementTypesVector = TypesVector<Element>;

// Template Specialisations
template<>
SpinTypesVector::TypesVector(unsigned long numberOfAlphaElectrons,
                             unsigned long numberOfBetaElectrons);

template<>
SpinTypesVector::TypesVector(std::vector<Spin> types);

template<>
ElementTypesVector::TypesVector(std::vector<Element> types);



namespace YAML {
    template<typename Type> struct convert<TypesVector<Type>> {
    static Node encode(const TypesVector<Type> & tv){
        Node node;
        for (unsigned long i = 0; i < tv.numberOfEntities(); i++) {
            node.push_back(tv[i]);
        }
        return node;
    }

    static bool decode(const Node& nodes, TypesVector<Type> & rhs){
        if(!nodes.IsSequence())
            return false;

        TypesVector<Type> tv;
        for (const auto &i : nodes)
            tv.append(i.as<Type>());

        rhs = tv;
        return true;
    }
};

template<typename Type>
Emitter& operator<< (Emitter& out, const TypesVector<Type>& tv){
    out << YAML::Flow << BeginSeq;
    for (unsigned i = 0; i < tv.numberOfEntities(); ++i)
        out << tv[i];
    out << EndSeq;
    return out;
};
}



#endif //AMOLQCPP_TYPESVECTOR_H
