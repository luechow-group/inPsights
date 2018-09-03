//
// Created by Michael Heuer on 31.08.18.
//

#ifndef AMOLQCPP_ISLICEABLE_H
#define AMOLQCPP_ISLICEABLE_H

#include "ReturnAndReset.h"
#include "Interval.h"
#include <memory>
#include <Eigen/Core>
#include "AbstractVector.h"


template<typename Scalar>
class ISliceable : public AbstractVector {
protected:

    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic,1>;
    using RefVectorType = Eigen::Ref<VectorType>;

public:

    ISliceable(long numberOfEntities = 0, long entityLength = 1)
            :
            AbstractVector(numberOfEntities, entityLength),
            data_(VectorType::Constant(numberOfEntities, 0)),
            resetType_(Reset::Automatic),
            sliceInterval_({0, numberOfEntities}),
            refPtr_(makeRefPtr(sliceInterval_))
            {}

    ISliceable(const ISliceable &rhs)
            :
            AbstractVector(rhs.numberOfEntities(), rhs.entityLength()),
            data_(rhs.asEigenVector()),
            resetType_(Reset::Automatic),
            sliceInterval_({0, rhs.numberOfEntities()}),
            refPtr_(makeRefPtr(sliceInterval_))
            {}

    ISliceable &operator=(const ISliceable &rhs) {
        if (this == &rhs) {
            resetRef();
            return *this;
        }
        AbstractVector::setEntityLength(rhs.entityLength());
        AbstractVector::setNumberOfEntities(rhs.numberOfEntities());
        data_ = rhs.asEigenVector();
        resetRef();

        return *this;
    }

    std::unique_ptr<RefVectorType> makeAllRefPtr(){
        return std::make_unique<RefVectorType>(data_.segment(0,numberOfEntities()));
    }
    std::unique_ptr<RefVectorType> makeRefPtr(const Interval& interval){
        return std::make_unique<RefVectorType>(data_.segment(calculateIndex(interval.start()),interval.numberOfEntities()));
    }

    void resetRef() {
        resetType_ = Reset::Automatic;
        sliceInterval_ = {0, numberOfEntities()};
        refPtr_.reset();
        refPtr_ = makeAllRefPtr();
    }

    void slice(const Interval& interval, const Reset& resetType = Reset::Automatic) {
        assert(interval.numberOfEntities() <= numberOfEntities() && "The interval is too long.");
        resetType_ = resetType;
        sliceInterval_ = interval;
        refPtr_.reset();
        refPtr_ = makeRefPtr(interval);
    }

    void entity(long i, const Reset& resetType = Reset::Automatic) {
        return slice(Interval(i), resetType);
    }


    void resetStrategy(const Usage &usage) {
        if( resetType_ == Reset::Automatic
            || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
            resetRef();
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        permuteMethod(permutation);
        resetRef();
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const Usage &usage) {
        permuteMethod(permutation);
        resetStrategy(usage);
    }

    void permuteMethod(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
        assert(permutation.indices().size() == sliceInterval_.numberOfEntities()
               && "The permutation vector length must be equal to the number of entities");

        auto tmp = resetType_;
        resetType_ = Reset::OnFinished;

        dataRef(Usage::NotFinished) = permutation* dataRef(Usage::NotFinished);

        resetType_ = tmp;
    }

    const Eigen::VectorXi& asEigenVector() const {
        return data_;
    }

    Eigen::VectorXi& asEigenVector() {
        return data_;
    }

    RefVectorType dataRef(const Usage &usage = Usage::NotFinished){
        if( resetType_ == Reset::Automatic
            || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
            return RETURN_AND_RESET<ISliceable<Scalar>,RefVectorType>(*this,*refPtr_).returnAndReset();
        else
            return *refPtr_;
    }


    bool operator==(const ISliceable<Scalar> &other) const {
        return (data_ == other.data_) && (sliceInterval_ == other.getSliceInterval())
               && (resetType_ == other.getResetType());
    }

    bool operator!=(const ISliceable<Scalar> &other) const {
        return !(*this == other);
    }

protected:
    const Reset &getResetType() const { return resetType_; }
    const Interval &getSliceInterval() const { return sliceInterval_; }


    VectorType data_;
    Reset resetType_;
    Interval sliceInterval_;
    std::unique_ptr<RefVectorType> refPtr_;

};

#endif //AMOLQCPP_ISLICEABLE_H
