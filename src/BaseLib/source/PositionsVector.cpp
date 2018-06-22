//
// Created by Michael Heuer on 29.10.17.
//

#include "PositionsVector.h"
#include "ToString.h"
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

using namespace Eigen;

PositionsVector::PositionsVector()
        : AbstractVector(),
          positions_(0)
{}

PositionsVector::PositionsVector(const VectorXd &positions) {
    auto size = positions.size();
    assert(size >= 0 && "Vector cannot be empty");
    assert(size%3 == 0 && "Vector must be 3N-dimensional");

    AbstractVector::setNumberOfEntities(size/3);
    positions_ = positions;
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
}

std::ostream& operator<<(std::ostream& os, const PositionsVector& pc){
    for (unsigned long i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::unsignedLongToString(i + 1) << " " << ToString::vector3dToString(pc[i]) << std::endl;
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

void PositionsVector::permute(long i, long j) {
    if(i != j) {
        Eigen::Vector3d temp = positions_.segment(calculateIndex(i),entityLength_);
        positions_.segment(calculateIndex(i),entityLength_) = positions_.segment(calculateIndex(j),entityLength_);
        positions_.segment(calculateIndex(j),entityLength_) = temp;
    }
}

long PositionsVector::calculateIndex(long i) const {
    return AbstractVector::calculateIndex(i)*entityLength_;
}

Eigen::Ref<Eigen::Vector3d> PositionsVector::operator()(long i){
    return Eigen::Ref<Eigen::Vector3d>(positions_.segment(i*entityLength_,entityLength_));
}

const Eigen::Ref<const Eigen::Vector3d>& PositionsVector::operator()(long i) const{
    return Eigen::Ref<const Eigen::Vector3d>(positions_.segment(i*entityLength_,entityLength_));
}

namespace YAML {
    Node convert<PositionsVector>::encode(const PositionsVector &rhs) {
        Node node;

        for (unsigned i = 0; i < rhs.numberOfEntities(); ++i) {
            node.push_back(rhs[i]);
        }
        return node;
    }
    bool convert<PositionsVector>::decode(const Node &node, PositionsVector &rhs) {
        if (!node.IsScalar())
            return false;
        PositionsVector pv;
        for (unsigned i = 0; i < node.size(); ++i)
            pv[i] = node[i].as<Eigen::Vector3d>();
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