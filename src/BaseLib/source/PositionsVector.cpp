//
// Created by Michael Heuer on 29.10.17.
//

#include <PositionsVector.h>
#include <ToString.h>
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>
#include <PositionsVector.h>
#include <PositionsVectorTransformer.h>


using namespace Eigen;

PositionsVector::PositionsVector()
        : AbstractVector(),
          positions_(0),
          resetType_(Reset::Automatic),
          sliceInterval_({0,0}),
          positionsRefPtr_(std::make_unique<PositionsRef>(positions_))
{}

PositionsVector::PositionsVector(const VectorXd &positions)
: PositionsVector() {
    auto size = positions.size();
    assert(size >= 0 && "Vector cannot be empty");
    assert(size%3 == 0 && "Vector must be 3N-dimensional");

    AbstractVector::setNumberOfEntities(size/3);
    positions_ = positions;

    resetRef();
}

Eigen::Vector3d PositionsVector::operator[](long i) const {
    return positions_.segment(calculateIndex(i),entityLength_);
}

Eigen::Vector3d PositionsVector::position(long i, const Usage& usage) {
    if( resetType_ == Reset::Automatic
    || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
        return RETURN_AND_RESET<PositionsVector,Eigen::Vector3d>(*this,positionsRef().segment(calculateIndex(i),entityLength_)).returnAndReset();
    else
        return positionsRef().segment(calculateIndex(i),entityLength_);
}

void PositionsVector::insert(const Eigen::Vector3d &position, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    long start = i*entityLength_;

    VectorXd before = positions_.head(start);
    VectorXd after = positions_.tail(numberOfEntities()*entityLength_-start);

    positions_.resize(numberOfEntities()*entityLength_+entityLength_);
    positions_ << before, position, after;

    incrementNumberOfEntities();
    resetRef();
}

std::ostream& operator<<(std::ostream& os, const PositionsVector& pc){
    for (long i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::longToString(i + 1) << " " << ToString::vector3dToString(pc[i]) << std::endl;
    }
    const_cast<PositionsVector&>(pc).resetRef(); //TODO refactor
    return os;
}

void PositionsVector::prepend(const Eigen::Vector3d &position) {
    this->insert(position,0);
}

void PositionsVector::append(const Eigen::Vector3d &position) {
    this->insert(position,numberOfEntities());
}

const Eigen::VectorXd & PositionsVector::positionsAsEigenVector() const {
    return positions_;
}

Eigen::VectorXd & PositionsVector::positionsAsEigenVector() {
    return positions_;
}

PositionsRef PositionsVector::positionsRef(const Usage& usage){
    if( resetType_ == Reset::Automatic
        || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
        return RETURN_AND_RESET<PositionsVector,PositionsRef>(*this,*positionsRefPtr_).returnAndReset();
    else
    return *positionsRefPtr_;
}

void PositionsVector::resetRef() {
    resetType_ = Reset::Automatic;
    positionsRefPtr_.reset();
    positionsRefPtr_ = std::make_unique<PositionsRef>(positions_);
    sliceInterval_ = {0,numberOfEntities()};
}

void PositionsVector::resetStrategy(const Usage &usage) {
    if( resetType_ == Reset::Automatic
        || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
        resetRef();
}

Eigen::PermutationMatrix<Eigen::Dynamic> PositionsVector::adaptedToEntityLength(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation){
    Eigen::VectorXi raw(permutation.indices().size()*entityLength_);

    for (int i = 0; i < permutation.indices().size(); ++i) {
        auto originIdx = i*entityLength_;
        auto targetIdx = permutation.indices()[i]*entityLength_;

        raw[originIdx+0] = targetIdx+0;
        raw[originIdx+1] = targetIdx+1;
        raw[originIdx+2] = targetIdx+2;
    }
    return PermutationMatrix<Eigen::Dynamic>(raw);
}

void PositionsVector::permuteMethod(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
    assert(permutation.indices().size() == sliceInterval_.numberOfEntities()
           && "The permutation vector length must be equal to the number of entities");

    auto tmp = resetType_;
    resetType_=Reset::OnFinished;

    positionsRef(Usage::NotFinished) = adaptedToEntityLength(permutation)*positionsRef(Usage::NotFinished);

    resetType_ = tmp;
}

void PositionsVector::permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const Usage &usage) {
    permuteMethod(permutation);
    resetStrategy(usage);
}

void PositionsVector::permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
    permuteMethod(permutation);
    resetRef();
}

void PositionsVector::translate(const Eigen::Vector3d &shift, const Usage& usage) {
    auto tmp = resetType_;
    resetType_ = Reset::OnFinished;

    for (long i = 0; i < sliceInterval_.numberOfEntities(); ++i)
        positionsRef(Usage::NotFinished).segment(calculateIndex(i),3) += shift;

    resetType_ = tmp;
    resetStrategy(usage);
}

void PositionsVector::rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection, const Usage& usage) {
    rotate(angle,{0,0,0},axisDirection,usage);
}

void PositionsVector::rotate(double angle,const Eigen::Vector3d &center,const Eigen::Vector3d &axisDirection, const Usage& usage) {
    auto rotMat = PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle, axisDirection - center));

    auto tmp = resetType_;
    resetType_ = Reset::OnFinished;

    this->translate(-center,Usage::NotFinished);

    for (long i = 0; i < sliceInterval_.numberOfEntities(); ++i)
        positionsRef().segment(calculateIndex(i), entityLength_) = position(i,Usage::NotFinished).transpose() * rotMat;

    this->translate(center,Usage::NotFinished);

    resetType_ = tmp;
    resetStrategy(usage);
};

PositionsVector::PositionsVector(const PositionsVector& rhs)
        : AbstractVector(rhs) {
    positions_ = rhs.positionsAsEigenVector();
    resetRef();
}

PositionsVector& PositionsVector::operator=(const PositionsVector& rhs){
    if(this == &rhs) {
        resetRef();
        return *this;
    }

    AbstractVector::setNumberOfEntities(rhs.numberOfEntities());
    positions_ = rhs.positionsAsEigenVector();
    resetRef();

    return *this;
}

PositionsVector& PositionsVector::entity(long i, const Reset& resetType) {
    return slice(Interval(i), resetType);
}

PositionsVector& PositionsVector::slice(const Interval& interval, const Reset& resetType) {
    assert(interval.numberOfEntities() <= numberOfEntities() && "The interval is too long.");
    resetType_ = resetType;
    sliceInterval_ = interval;
    positionsRefPtr_.reset();
    positionsRefPtr_ = std::make_unique<PositionsRef>(
            positions_.segment(calculateIndex(interval.start()),interval.numberOfEntities()*entityLength_));
    return *this;
}

bool PositionsVector::operator==(const PositionsVector& other) const {
    return (positions_.isApprox(other.positions_),0)
    && (sliceInterval_ == other.getSliceInterval())
    && (resetType_ == other.getResetType());
}

bool PositionsVector::operator!=(const PositionsVector&other) const {
    return !(*this == other);
}

long PositionsVector::calculateIndex(long i) const {
    return AbstractVector::calculateIndex(i)*entityLength_;
}

namespace YAML {
    Node convert<PositionsVector>::encode(const PositionsVector &rhs) {
        Node node;
        for (unsigned i = 0; i < rhs.numberOfEntities(); ++i)
            node.push_back(rhs[i]);
        return node;
    }
    bool convert<PositionsVector>::decode(const Node &node, PositionsVector &rhs) {
        if (!node.IsSequence())
            return false;
        PositionsVector pv;
        for (const auto &i : node)
            pv.append(i.as<Eigen::Vector3d>());
        rhs = pv;
        return true;
    }

    Emitter &operator<<(Emitter &out, const PositionsVector &p) {
        out << Flow << BeginSeq;
        for (unsigned i = 0; i < p.numberOfEntities(); ++i)
            out << p[i];
        out << EndSeq;
        return out;
    }
}
