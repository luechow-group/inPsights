// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <PositionsVector.h>
#include <ToString.h>
#include <EigenYamlConversion.h>
#include <PositionsVectorTransformer.h>
#include <yaml-cpp/yaml.h>
#include <random>

using namespace Eigen;

PositionsVector::PositionsVector()
        : InsertableVector(0,3) {}

PositionsVector::PositionsVector(const VectorXd &positions)
: PositionsVector() {
    auto size = positions.size();
    assert(size >= 0 && "Vector cannot be empty");
    assert(size%entityLength() == 0 && "Vector must be 3N-dimensional");

    AbstractVector::setNumberOfEntities(size/entityLength());
    data_ = positions;
}

Eigen::Vector3d PositionsVector::operator[](long i) const {
    return data_.segment(calculateIndex(i),entityLength());
}

Eigen::Vector3d PositionsVector::position(long i) {
    return operator[](i);
}

void PositionsVector::prepend(const Eigen::Vector3d &position) {
    insert(position,0);
}
void PositionsVector::append(const Eigen::Vector3d &position) {
    insert(position,AbstractVector::numberOfEntities());
}
void PositionsVector::insert(const Eigen::Vector3d &position, long i) {
    InsertableVector<double>::insert(position,i);
}

void PositionsVector::translate(const Eigen::Vector3d &shift) {
   for (long i = 0; i < numberOfEntities(); ++i)
        data_.segment(calculateIndex(i), 3) += shift;
}

void PositionsVector::rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection) {
    rotate(angle,{0,0,0},axisDirection);
}

void PositionsVector::rotate(double angle,const Eigen::Vector3d &center,const Eigen::Vector3d &axisDirection) {
    auto rotMat = PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle, axisDirection));

    this->translate(-center);

    for (long i = 0; i < numberOfEntities(); ++i)
        data_.segment(calculateIndex(i), entityLength()) = position(i).transpose() * rotMat;

    this->translate(center);
};

std::ostream& operator<<(std::ostream& os, const PositionsVector& pc){
    for (Eigen::Index i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::longToString(i + 1) << " " << ToString::vector3dToString(pc[i]) << std::endl;
    }
    return os;
}

bool PositionsVector::operator==(const PositionsVector& other) const {
    return DataVector<double>::operator==(other);
}

bool PositionsVector::operator!=(const PositionsVector&other) const {
    return !(*this == other);
}


void PositionsVector::shake(double radius, std::default_random_engine& rng){
    auto maxDev = 1./std::sqrt(3.0)*radius;
    std::uniform_real_distribution<double> uniformRealDistribution(-maxDev, maxDev);

    auto eLength = entityLength();
    for (long i = 0; i < numberOfEntities(); ++i) {
        data_[i*eLength+0] += uniformRealDistribution(rng);
        data_[i*eLength+1] += uniformRealDistribution(rng);
        data_[i*eLength+2] += uniformRealDistribution(rng);
    }
};


namespace YAML {
    Node convert<PositionsVector>::encode(const PositionsVector &rhs) {
        Node node;
        for (long i = 0; i < rhs.numberOfEntities(); ++i)
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
        out << Flow << BeginSeq << Newline;
        for (long i = 0; i < p.numberOfEntities(); ++i)
            out << p[i] << Newline;
        out << EndSeq << Auto;
        return out;
    }
}
