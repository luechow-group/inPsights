//
// Created by Michael Heuer on 22.04.18.
//

#ifndef INPSIGHTS_TYPESVECTOR_H
#define INPSIGHTS_TYPESVECTOR_H

#include "SpinType.h"
#include "ElementType.h"
#include "EnumeratedType.h"
#include "InsertableVector.h"
#include <vector>
#include <yaml-cpp/yaml.h>

template <typename Type>
class TypesVector : public InsertableVector<int> {
public:

    // unsigned long dummy is there to be specialized in TypesVector<Spin::SpinTypes>
    TypesVector(unsigned long size = 0, unsigned long dummy = 0)
    : InsertableVector<int>(long(size)) {}

    TypesVector(const Eigen::VectorXi &types)
            : TypesVector() {
        auto size = types.size();
        assert(size >= 0 && "Vector cannot be empty");

        AbstractVector::setNumberOfEntities(size);
        data_ = types;
    }


    TypesVector(std::vector<Type> types)
    : InsertableVector<int>(0)
    {
        for (const auto& type : types){
            assert(int(type) >= int(Spins::first()));
            assert(int(type) <= int(Elements::last()));
            this->append(type);
        }
    }

    void prepend(const Type& type) { insert(type,0); }
    void append(const Type& type) { insert(type,numberOfEntities()); }
    void insert(const Type& type, long i) {
        Eigen::Matrix<int,1,1>elem;
        elem << int(type);
        InsertableVector<int>::insert(elem,i);
    }

    Type type(long i) {
        return operator[](i);
    }

    Type operator[](long i) const {
        return static_cast<Type>(data_[calculateIndex(i)]);
    }
    
    bool operator==(const TypesVector<Type> &other) const {
        return DataVector<int>::operator==(other);
    }

    bool operator!=(const TypesVector<Type> &other) const {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const TypesVector& tv){
        for (long i = 0; i < tv.numberOfEntities(); i++) {
            os << tv[i] << std::endl;
        }
        return os;
    }

    unsigned countOccurence(const Type &type) const { // cannot handle slices
        unsigned count = 0;
        for (long i = 0; i < this->numberOfEntities(); ++i) {
            if (this->operator[](i) == type)
                count++;
        }
        return count;
    }

    EnumeratedType<Type> getEnumeratedTypeByIndex(long i) const { // cannot handle slices
        auto type = this->operator[](i);
        unsigned count = 0;

        for (long j = 0; j < calculateIndex(i); ++j)
            if(this->operator[](j) == type)
                count++;

        return {type,count};
    };

    std::pair<bool,long> findIndexOfEnumeratedType(EnumeratedType<Type> indexedType) const { // cannot handle slices
        assert(indexedType.number_ < numberOfEntities() &&
               "This index is out of bounds.");

        unsigned count = 0;
        for (long j = 0; j < numberOfEntities(); ++j) {
            if(this->operator[](j) == indexedType.type_) count++;
            if(count-1 == indexedType.number_) return {true, j};
        }
        return {false,0};
    }

    //TODO dirty - refactor
    std::vector<EnumeratedType<Type>> countTypes() const {

        std::vector<EnumeratedType<Type>> typeCountsPair;
        if(numberOfEntities() > 0) {
            if(numberOfEntities() == 1) {
                return {{this->operator[](0),1}};
            }
            auto copy = data_;
            std::sort(copy.data(), copy.data()+numberOfEntities());

            unsigned count = 1;
            for (long i = 1; i < numberOfEntities(); ++i) {
                const auto& lastType = copy[i - 1];
                const auto& currentType = copy[i];

                if (currentType == lastType)
                    count++;
                else {
                    typeCountsPair.emplace_back(EnumeratedType<Type>(Type(lastType), count));
                    count = 1;
                }
                // treat last element
                if(i == numberOfEntities()-1)
                    typeCountsPair.emplace_back(EnumeratedType<Type>(Type(currentType), count));
            }
        }
        return typeCountsPair;
    }

    unsigned multiplicity() = delete;
    void flipSpins() = delete;
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
unsigned SpinTypesVector::multiplicity();

template<>
void SpinTypesVector::flipSpins();

template<>
ElementTypesVector::TypesVector(std::vector<Element> types);

namespace YAML {
    template<typename Type> struct convert<TypesVector<Type>> {
    static Node encode(const TypesVector<Type> & tv){
        Node node;
        for (long i = 0; i < tv.numberOfEntities(); i++) {
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



#endif //INPSIGHTS_TYPESVECTOR_H
