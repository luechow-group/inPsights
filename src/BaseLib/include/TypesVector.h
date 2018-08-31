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

using TypesRef = Eigen::Ref<Eigen::VectorXi>;

template <typename Type>
class TypesVector : public AbstractVector {
public:

    // unsigned long empty is there to be specialized in TypesVector<Spin::SpinTypes>
    TypesVector(unsigned long size = 0, unsigned long empty = 0)
            : AbstractVector(size),
              types_(Eigen::VectorXi::Constant(size, 0)),
              resetType_(Reset::Automatic),
              sliceInterval_({0,numberOfEntities()}),
              typesRefPtr_(std::make_unique<TypesRef>(types_))
              {}

    TypesVector(std::vector<Type> types)
            : AbstractVector(0),
              types_(0),
              resetType_(Reset::Automatic),
              sliceInterval_({0,0}),
              typesRefPtr_()
    {
        for (const auto& type : types){
            assert(int(type) >= int(Spins::first()));
            assert(int(type) <= int(Elements::last()));
            this->append(type);
        }
        resetRef();
    }

    TypesVector(const TypesVector& rhs)
    : AbstractVector(rhs) {
        types_ = rhs.typesAsEigenVector();
        resetRef();
    }

    TypesVector& operator=(const TypesVector& rhs){
        if(this == &rhs) {
            resetRef();
            return *this;
        }

        AbstractVector::setNumberOfEntities(rhs.numberOfEntities());
        types_ = rhs.typesAsEigenVector();
        resetRef();

        return *this;
    }


    TypesVector<Type>& slice(const Interval& interval, const Reset& resetType = Reset::Automatic) {
        assert(interval.numberOfEntities() <= numberOfEntities() && "The interval is too long.");
        resetType_ = resetType;
        sliceInterval_ = interval;
        typesRefPtr_.reset();
        typesRefPtr_ = std::make_unique<TypesRef>(
                types_.segment(calculateIndex(interval.start()),interval.numberOfEntities()));
        return *this;
    }

    TypesVector<Type>& entity(long i, const Reset& resetType = Reset::Automatic) {
        return slice(Interval(i), resetType);
    }


    TypesRef typesRef(const Usage& usage = Usage::NotFinished){
        if( resetType_ == Reset::Automatic
            || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
            return RETURN_AND_RESET<TypesVector<Type>,TypesRef>(*this,*typesRefPtr_).returnAndReset();
        else
            return *typesRefPtr_;
    }

    void resetRef() {
        resetType_ = Reset::Automatic;
        typesRefPtr_.reset();
        typesRefPtr_ = std::make_unique<TypesRef>(types_);
        sliceInterval_ = {0,numberOfEntities()};
    }

    void resetStrategy(const Usage &usage) {
        if( resetType_ == Reset::Automatic
            || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
            resetRef();
    }

    Type operator[](long i) const {
        return static_cast<Type>(types_[calculateIndex(i)]);
    }

    void insert(Type type, long i) {
        assert(i >= 0 && "The index must be positive.");
        assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

        Eigen::VectorXi before = types_.head(i);
        Eigen::VectorXi after = types_.tail(numberOfEntities()-i);

        types_.resize(numberOfEntities()+1);
        types_ << before, int(type), after;

        incrementNumberOfEntities();
        resetRef();
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

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        permuteMethod(permutation);
        resetRef();
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const Usage &usage) {
        permuteMethod(permutation);
        resetStrategy(usage);
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
            auto copy = types_;
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

    const Reset& getResetType() const { return resetType_; }
    const Interval& getSliceInterval() const { return sliceInterval_; }

    bool operator==(const TypesVector<Type>& other) const {
        return (types_ == other.types_)
        && (sliceInterval_ == other.getSliceInterval())
        && (resetType_ == other.getResetType());
    }

    bool operator!=(const TypesVector<Type>& other) const {
        return !(*this == other);
    }

private:
    Eigen::VectorXi types_;

    Reset resetType_;
    Interval sliceInterval_;
    std::unique_ptr<TypesRef> typesRefPtr_;

    void permuteMethod(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
        assert(permutation.indices().size() == sliceInterval_.numberOfEntities()
               && "The permutation vector length must be equal to the number of entities");

        auto tmp = resetType_;
        resetType_ = Reset::OnFinished;

        typesRef(Usage::NotFinished) = permutation*typesRef(Usage::NotFinished);

        resetType_ = tmp;
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
