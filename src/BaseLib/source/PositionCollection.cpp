//
// Created by Michael Heuer on 29.10.17.
//

#include "PositionCollection.h"
#include "ToString.h"

using namespace Eigen;

PositionCollection::PositionCollection()
        : AbstractCollection(),
          positions_(0)
{}

PositionCollection::PositionCollection(const VectorXd &positions) {
    auto size = positions.size();
    assert(size >= 0 && "Vector cannot be empty");
    assert(size%3 == 0 && "Vector must be 3N-dimensional");

    AbstractCollection::setNumberOfEntities(size/3);
    positions_ = positions;
}

Eigen::Vector3d PositionCollection::operator[](long i) const {
    return positions_.segment(calculateIndex(i),entityLength_);
}

void PositionCollection::insert(const Eigen::Vector3d &position, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    long start = i*entityLength_;

    VectorXd before = positions_.head(start);
    VectorXd after = positions_.tail(numberOfEntities()*entityLength_-start);

    positions_.resize(numberOfEntities()*entityLength_+entityLength_);
    positions_ << before, position, after;

    incrementNumberOfEntities();
}

std::ostream& operator<<(std::ostream& os, const PositionCollection& pc){
    for (unsigned long i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::vector3d2string(pc[i]) << std::endl;
    }
    return os;
}

void PositionCollection::prepend(const Eigen::Vector3d &position) {
    this->insert(position,0);
}

void PositionCollection::append(const Eigen::Vector3d &position) {
    this->insert(position,numberOfEntities());
}

const Eigen::VectorXd & PositionCollection::positionsAsEigenVector() const {
    return positions_;
}

Eigen::VectorXd & PositionCollection::positionsAsEigenVector() {
    return positions_;
}

void PositionCollection::permute(long i, long j) {
    if(i != j) {
        Eigen::Vector3d temp = positions_.segment(calculateIndex(i),entityLength_);
        positions_.segment(calculateIndex(i),entityLength_) = positions_.segment(calculateIndex(j),entityLength_);
        positions_.segment(calculateIndex(j),entityLength_) = temp;
    }
}

long PositionCollection::calculateIndex(long i) const {
    return AbstractCollection::calculateIndex(i)*entityLength_;
}
