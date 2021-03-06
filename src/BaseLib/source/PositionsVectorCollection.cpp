// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <PositionsVectorCollection.h>
#include <yaml-cpp/yaml.h>


PositionsVectorCollection::PositionsVectorCollection()
        : AbstractVector(0),
          positionsVectorCollection_(0),
          numberOfPositionEntities_(0)
{}

PositionsVectorCollection::PositionsVectorCollection(const std::vector<PositionsVector> &positionsVectorCollection)
        : AbstractVector(positionsVectorCollection.size()),
          numberOfPositionEntities_(0)
{
    if( !positionsVectorCollection.empty()) {
        numberOfPositionEntities_ = positionsVectorCollection[0].numberOfEntities();

        assert(std::all_of(std::next(positionsVectorCollection.cbegin()), positionsVectorCollection.cend(),
                [this](const PositionsVector &p){ return p.numberOfEntities() == numberOfPositionEntities_; })
                && "All position collections must contain the same number of positions.");
    }
}

long PositionsVectorCollection::numberOfPositionsEntities() const {
    return numberOfPositionEntities_;
}

PositionsVector PositionsVectorCollection::operator[](long i) const {
    return positionsVectorCollection_[calculateIndex(i)];
}

void PositionsVectorCollection::insert(const PositionsVector &positionsVector, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");
    if(numberOfEntities() == 0) {
        numberOfPositionEntities_ = positionsVector.numberOfEntities();
    }
    else {
        assert(numberOfPositionEntities_ == positionsVector.numberOfEntities()
               && "All position collections must contain the same number of positions.");
    }

    positionsVectorCollection_.insert(positionsVectorCollection_.begin()+i,positionsVector);
    incrementNumberOfEntities();
}

void PositionsVectorCollection::append(const PositionsVector &positionsVector) {
    insert(positionsVector,numberOfEntities());
}

void PositionsVectorCollection::prepend(const PositionsVector &positionsVector) {
    insert(positionsVector,0);
}

void PositionsVectorCollection::permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
    for(auto & pv : positionsVectorCollection_){
        pv.permute(permutation);
    }
}

double PositionsVectorCollection::norm(long i, long j) const{
    return (positionsVectorCollection_[i].asEigenVector()
            - positionsVectorCollection_[j].asEigenVector()).norm();
}

const std::vector<PositionsVector>& PositionsVectorCollection::positionsVectorCollection() const {
    return positionsVectorCollection_;
}

std::vector<PositionsVector>& PositionsVectorCollection::positionsVectorCollection() {
    return positionsVectorCollection_;
}


namespace YAML {
    Node convert<PositionsVectorCollection>::encode(const PositionsVectorCollection &rhs) {
        Node node;
        for (unsigned i = 0; i < rhs.numberOfEntities(); ++i)
            node.push_back(rhs[i]);
        return node;
    }
    bool convert<PositionsVectorCollection>::decode(const Node &node, PositionsVectorCollection &rhs) {
        if (!node.IsScalar())
            return false;
        PositionsVectorCollection pvc;
        for (unsigned i = 0; i < node.size(); ++i)
            pvc[i] = node[i].as<PositionsVector>();
        rhs = pvc;
        return true;
    }

    Emitter &operator<<(Emitter &out, const PositionsVectorCollection &pvc) {
        out << Flow << BeginSeq << Newline;
        for (unsigned i = 0; i < pvc.numberOfEntities(); ++i)
            out << pvc[i] << Newline;
        out << EndSeq;
        return out;
    }
}