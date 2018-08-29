//
// Created by Michael Heuer on 29.10.17.
//

#include "PositionsVector.h"
#include "ToString.h"
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>
#include <PositionsVector.h>
#include <PositionsVectorTransformer.h>

using namespace Eigen;

PositionsVector::PositionsVector()
        : AbstractVector(),
          positions_(0),
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

    resetRefToAll();
}

Eigen::Vector3d PositionsVector::operator[](long i) const {
    return positions_.segment(calculateIndex(i),entityLength_);
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
    resetRefToAll();
}

std::ostream& operator<<(std::ostream& os, const PositionsVector& pc){
    for (long i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::longToString(i + 1) << " " << ToString::vector3dToString(pc[i]) << std::endl;
    }
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

PositionsRef PositionsVector::positionsRef(){
    return *positionsRefPtr_;
}

void PositionsVector::resetRefToAll() {
    positionsRefPtr_.reset();
    positionsRefPtr_ = std::make_unique<PositionsRef>(positions_);
    sliceInterval_ = {0,numberOfEntities()};
}

void PositionsVector::permute(long i, long j) {
    if(i != j) {
        Eigen::Vector3d temp = positions_.segment(calculateIndex(i),entityLength_);
        positions_.segment(calculateIndex(i),entityLength_) = positions_.segment(calculateIndex(j),entityLength_);
        positions_.segment(calculateIndex(j),entityLength_) = temp;
    }
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

void PositionsVector::permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
    assert(permutation.indices().size() == sliceInterval_.numberOfEntities()
    && "The permutation vector length must be equal to the number of entities");

    positionsRef() = adaptedToEntityLength(permutation)*positionsRef();

    resetRefToAll();
}

void PositionsVector::translate(const Eigen::Vector3d &shift, bool resetAfterwardsQ) {

    for (long i = 0; i < sliceInterval_.numberOfEntities(); ++i)
        positionsRef().segment(calculateIndex(i),3) += shift;

    if(resetAfterwardsQ) resetRefToAll();
}

void PositionsVector::rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection, bool resetAfterwardsQ) {
    rotate(angle,{0,0,0},axisDirection,resetAfterwardsQ);
}

void PositionsVector::rotate(double angle,const Eigen::Vector3d &center,const Eigen::Vector3d &axisDirection, bool resetAfterwardsQ) {

    auto rotMat = PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle, axisDirection - center));

    this->translate(-center,false);

    for (long i = 0; i < sliceInterval_.numberOfEntities(); ++i) {
        long idx = sliceInterval_.start()+i;

        //TODO refactoring needed: implement methods to get the i-th Vector3d of a slice
        positionsRef().segment(calculateIndex(i), entityLength_) = this->operator[](idx).transpose() * rotMat;
    }
    this->translate(center,false);

    if (resetAfterwardsQ) resetRefToAll();
};

PositionsVector::PositionsVector(const PositionsVector& rhs)
        : AbstractVector(rhs) {
    positions_ = rhs.positionsAsEigenVector();
    resetRefToAll();
}

PositionsVector& PositionsVector::operator=(const PositionsVector& rhs){
    if(this == &rhs) {
        resetRefToAll();
        return *this;
    }

    AbstractVector::setNumberOfEntities(rhs.numberOfEntities());
    positions_ = rhs.positionsAsEigenVector();
    resetRefToAll();

    return *this;
}

PositionsVector& PositionsVector::entity(long i) {
    sliceInterval_ = Interval(i);
    return slice(sliceInterval_);
}

PositionsVector& PositionsVector::slice(const Interval& interval) {
    assert(interval.checkBounds(numberOfEntities()) && "The end variable of the interval is out of bounds.");
    sliceInterval_ = interval;
    positionsRefPtr_.reset();
    positionsRefPtr_ = std::make_unique<PositionsRef>(
            positions_.segment(calculateIndex(interval.start()),interval.numberOfEntities()*entityLength_));
    return *this;
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
