//
// Created by Michael Heuer on 22.04.18.
//

#ifndef AMOLQCPP_TYPESVECTOR_H
#define AMOLQCPP_TYPESVECTOR_H

#include <Eigen/Core>
#include "AbstractVector.h"
#include "SpinType.h"
#include "ElementType.h"
#include <vector>

template <typename Type>
class TypesVector : public AbstractVector {
public:

    // unsigned long empty is there to be specialized in TypesVector<Spin::SpinTypes>
    TypesVector(unsigned long size = 0, unsigned long empty = 0)
            : AbstractVector(size),
              types_(Eigen::VectorXi::Constant(size, 0))
    {}

    //TODO specialise
    TypesVector(std::vector<Type> types)
            : AbstractVector(0),
              types_(0)
    {
        for (const auto& type : types){
            assert(int(type) >= int(Spins::first()));
            assert(int(type) <= int(Elements::last()));
            this->append(type);
        }
    }

    /*TypesVector(const Eigen::VectorXi& typesAsEigenVector)
            : AbstractVector(typesAsEigenVector.size()),
              types_(typesAsEigenVector)
    {
        assert(types_.minCoeff() >= int(Spins::first()));
        assert(types_.maxCoeff() <= int(Elements::last()));
    }*/

    Type operator[](long i) const {
        return static_cast<Type>(types_[calculateIndex(i)]);
    }

    void insert(Type type, long i) {
        assert(i >= 0 && "The index must be positive.");
        assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

        //TODO! check if types are compatible
        //if ()

        Eigen::VectorXi before = types_.head(i);
        Eigen::VectorXi after = types_.tail(numberOfEntities()-i);

        types_.resize(numberOfEntities()+1);
        types_ << before, int(type), after;

        incrementNumberOfEntities();
    }

    void prepend(Type type) {
        this->insert(type,0);
    }

    void append(Type type) {
        this->insert(type,numberOfEntities());
    }

    const Eigen::VectorXi& typesAsEigenVector() const {
        return types_;
    }

    Eigen::VectorXi& typesAsEigenVector() {
        return types_;
    }

    void permute(long i, long j) {
        if(i != j) {
            int temp = types_[calculateIndex(i)];
            types_[calculateIndex(i)] = types_[calculateIndex(j)];
            types_[calculateIndex(j)] = temp;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const TypesVector& tv){
        for (unsigned long i = 0; i < tv.numberOfEntities(); i++) {
            os << tv[i] << std::endl;
        }
        return os;
    }

protected:
    Eigen::VectorXi types_;
};

using SpinTypesVector = TypesVector<Spins::SpinType>;
using ElementTypesVector = TypesVector<Elements::ElementType>;

// Template Specialisation
template<>
SpinTypesVector::TypesVector(unsigned long numberOfAlphaElectrons,
                             unsigned long numberOfBetaElectrons);

template<>
SpinTypesVector::TypesVector(std::vector<Spins::SpinType> types);

template<>
ElementTypesVector ::TypesVector(std::vector<Elements::ElementType> types);

#endif //AMOLQCPP_TYPESVECTOR_H
