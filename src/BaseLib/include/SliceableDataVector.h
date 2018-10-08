//
// Created by Michael Heuer on 31.08.18.
//

#ifndef AMOLQCPP_SLICEABLEDATAVECTOR_H
#define AMOLQCPP_SLICEABLEDATAVECTOR_H

#include "AbstractVector.h"
#include "ReturnAndReset.h"
#include "Interval.h"
#include <memory>
#include <Eigen/Core>

class ISliceable : public AbstractVector{
    friend class AbstractState;
public:
    void resetSlice(){
        resetType_ = Reset::Automatic;
        sliceInterval_ = {0, numberOfEntities()};
    }
protected:
    ISliceable(long numberOfEntities = 0, long entityLength = 1)
    :
    AbstractVector(numberOfEntities, entityLength),
    resetType_(Reset::Automatic),
    sliceInterval_({0, numberOfEntities})
    {}

    void setSlice(const Interval& interval, const Reset& resetType) {
        assert(interval.numberOfEntities() <= numberOfEntities() && "The interval is too long.");
        resetType_ = resetType;
        sliceInterval_ = interval;
    }

    bool operator==(const ISliceable& other) const{
        return (resetType_ == other.resetType())
        && (sliceInterval_ == other.sliceInterval());
    }

    bool operator!=(const ISliceable& other) const{
        return !(*this == other);
    }

    const Reset &resetType() const { return resetType_; }
    const Interval &sliceInterval() const { return sliceInterval_; }

Reset resetType_;
Interval sliceInterval_;
};

template<typename Scalar>
class SliceableDataVector : public ISliceable {

protected:
    using EigenVecType = Eigen::Matrix<Scalar, Eigen::Dynamic,1>;
    using RefEigenVecType = Eigen::Ref<EigenVecType>;

public:

    SliceableDataVector(long numberOfEntities = 0, long entityLength = 1)
    :
    ISliceable(numberOfEntities,entityLength),
    data_(EigenVecType::Constant(numberOfEntities, 0)),
    refPtr_()
    {
        resetRef();
    }

    SliceableDataVector(const SliceableDataVector &rhs)
    :
    ISliceable(rhs.numberOfEntities(), rhs.entityLength()),
    data_(rhs.asEigenVector()),
    refPtr_()
    {
        resetRef();
    }

    SliceableDataVector &operator=(const SliceableDataVector &rhs) {
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

    void assertResetInstruction(const Usage& usage){
        assert(
                ((resetType() == Reset::Automatic) && (usage == Usage::Standard)) ||
                ((resetType() == Reset::Manual) && (usage == Usage::Standard)) ||
                ((resetType() == Reset::OnFinished) && (usage == Usage::Finished || usage == Usage::NotFinished))
                && "You used conflicting reset instruction.");
    }

    void resetRef() {
        resetSlice();
        refPtr_.reset();
        refPtr_ = makeAllRefPtr();
    }

    void slice(const Interval& interval, const Reset& resetType = Reset::Automatic) {
        setSlice(interval,resetType);
        refPtr_.reset();
        refPtr_ = makeRefPtr(interval);
    }

    void entity(long i, const Reset& resetType = Reset::Automatic) {
        slice(Interval(i), resetType);
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        permuteMethod(permutation);
        resetRef();
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const Usage &usage) {
        permuteMethod(permutation);
        resetStrategy(usage);
    }

    const EigenVecType& asEigenVector() const {
        return data_;
    }

    EigenVecType& asEigenVector() {
        return data_;
    }

    RefEigenVecType dataRef(const Usage &usage = Usage::Standard){
        assertResetInstruction(usage);

        if((resetType() == Reset::Automatic && usage == Usage::Standard)
        || (resetType() == Reset::OnFinished && usage == Usage::Finished))
            return RETURN_AND_RESET<SliceableDataVector<Scalar>,RefEigenVecType>(*this,*refPtr_).returnAndReset();
        else
            return *refPtr_;
    }

protected:
    //TODO make double template?
    bool operator==(const SliceableDataVector<Scalar> &other) const {
        return (data_ == other.data_) && ISliceable::operator==(other);
    }

    //TODO make double template?
    bool operator!=(const SliceableDataVector<Scalar> &other) const {
        return !(*this == other);
    }

    void resetStrategy(const Usage &usage) {
        assertResetInstruction(usage);

        if((resetType() == Reset::Automatic && usage == Usage::Standard)
        || (resetType() == Reset::OnFinished && usage == Usage::Finished))
            resetRef();
    }

    std::unique_ptr<RefEigenVecType> makeAllRefPtr(){
        return std::make_unique<RefEigenVecType>(data_.segment(0,numberOfEntities()*entityLength()));
    }

    std::unique_ptr<RefEigenVecType> makeRefPtr(const Interval& interval){
        return std::make_unique<RefEigenVecType>(
                data_.segment(
                        calculateIndex(interval.start()),
                        interval.numberOfEntities()*entityLength())
        );
    }

    void permuteMethod(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
        assert(permutation.indices().size() == sliceInterval_.numberOfEntities()
               && "The permutation vector length must be equal to the number of entities");

        auto tmp = resetType();
        resetType_ = Reset::OnFinished; //TODO ADD SETTER

        if(entityLength() > 1) {
            dataRef(Usage::NotFinished) =  adaptedToEntityLength(permutation) * dataRef(Usage::NotFinished);
        }
        else
            dataRef(Usage::NotFinished) = permutation* dataRef(Usage::NotFinished);

        resetType_ = tmp; //TODO ADD SETTER
    }

    Eigen::PermutationMatrix<Eigen::Dynamic> adaptedToEntityLength(
            const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation){

        Eigen::VectorXi raw(permutation.indices().size()*entityLength());

        for (int i = 0; i < permutation.indices().size(); ++i) {
            auto originIdx = i*entityLength();
            auto targetIdx = permutation.indices()[i]*entityLength();

            raw[originIdx+0] = targetIdx+0;
            raw[originIdx+1] = targetIdx+1;
            raw[originIdx+2] = targetIdx+2;
        }
        return Eigen::PermutationMatrix<Eigen::Dynamic>(raw);
    }

    EigenVecType data_;
    std::unique_ptr<RefEigenVecType> refPtr_;
};

#endif //AMOLQCPP_SLICEABLEDATAVECTOR_H