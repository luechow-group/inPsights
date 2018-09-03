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
        : IInsertable(0,3)
{
    resetRef();
}

PositionsVector::PositionsVector(const VectorXd &positions)
: PositionsVector() {
    auto size = positions.size();
    assert(size >= 0 && "Vector cannot be empty");
    assert(size%3 == 0 && "Vector must be 3N-dimensional");

    AbstractVector::setNumberOfEntities(size/3);
    data_ = positions;

    resetRef();
}

Eigen::Vector3d PositionsVector::operator[](long i) const {
    return data_.segment(calculateIndex(i),entityLength());
}

// TODO REPLACE?
Eigen::Vector3d PositionsVector::position(long i, const Usage& usage) {
    if( resetType_ == Reset::Automatic
    || (resetType_ == Reset::OnFinished && usage == Usage::Finished))
        return RETURN_AND_RESET<PositionsVector,Eigen::Vector3d>(*this,dataRef().segment(calculateIndex(i),entityLength())).returnAndReset();
    else
        return dataRef().segment(calculateIndex(i),entityLength());
}


void PositionsVector::prepend(const Eigen::Vector3d &position) {
    insert(position,0);
}
void PositionsVector::append(const Eigen::Vector3d &position) {
    insert(position,AbstractVector::numberOfEntities());
}
void PositionsVector::insert(const Eigen::Vector3d &position, long i) {
    IInsertable<double>::insert(position,i);
}

void PositionsVector::translate(const Eigen::Vector3d &shift, const Usage& usage) {
    auto tmp = resetType_;
    resetType_ = Reset::OnFinished;

    for (long i = 0; i < sliceInterval_.numberOfEntities(); ++i)
        dataRef(Usage::NotFinished).segment(calculateIndex(i),3) += shift;

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
        dataRef().segment(calculateIndex(i), entityLength()) = position(i,Usage::NotFinished).transpose() * rotMat;

    this->translate(center,Usage::NotFinished);

    resetType_ = tmp;
    resetStrategy(usage);
};

PositionsVector& PositionsVector::entity(long i, const Reset& resetType) {
    return slice(Interval(i), resetType);
}

PositionsVector& PositionsVector::slice(const Interval& interval, const Reset& resetType) {
    ISliceable<double>::slice(interval,resetType);
    return *this;
}

std::ostream& operator<<(std::ostream& os, const PositionsVector& pc){
    for (long i = 0; i < pc.numberOfEntities(); i++){
        os << ToString::longToString(i + 1) << " " << ToString::vector3dToString(pc[i]) << std::endl;
    }
    const_cast<PositionsVector &>(pc).resetRef(); //TODO refactor
    return os;
}

bool PositionsVector::operator==(const PositionsVector& other) const {
    return ISliceable<double>::operator==(other);
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
        out << Flow << BeginSeq;
        for (unsigned i = 0; i < p.numberOfEntities(); ++i)
            out << p[i];
        out << EndSeq;
        return out;
    }
}
