//
// Created by Michael Heuer on 29.10.17.
//

#include <PositionsVector.h>
#include <ToString.h>
#include <EigenYamlConversion.h>
#include <PositionsVectorTransformer.h>
#include <yaml-cpp/yaml.h>

using namespace Eigen;

PositionsVector::PositionsVector()
        : InsertableVector(0,3) {}

PositionsVector::PositionsVector(const VectorXd &positions)
: PositionsVector() {
    auto size = positions.size();
    assert(size >= 0 && "Vector cannot be empty");
    assert(size%3 == 0 && "Vector must be 3N-dimensional");

    AbstractVector::setNumberOfEntities(size/3);
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
    for (long i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::longToString(i + 1) << " " << ToString::vector3dToString(pc[i]) << std::endl;
    }
    return os;
}

bool PositionsVector::operator==(const PositionsVector& other) const {
    return SliceableDataVector<double>::operator==(other);
}

bool PositionsVector::operator!=(const PositionsVector&other) const {
    return !(*this == other);
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
        out << Flow << BeginSeq << Newline;
        for (unsigned i = 0; i < p.numberOfEntities(); ++i)
            out << p[i] << Newline;
        out << EndSeq << Auto;
        return out;
    }
}
