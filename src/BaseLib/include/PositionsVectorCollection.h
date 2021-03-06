// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_POSITIONSVECTORCOLLECTION_H
#define INPSIGHTS_POSITIONSVECTORCOLLECTION_H

#include <vector>
#include "PositionsVector.h"

class PositionsVectorCollection : public AbstractVector{
public:
    PositionsVectorCollection();
    explicit PositionsVectorCollection(const std::vector<PositionsVector> &positionsVectorCollection);

    PositionsVector operator[](long i) const;

    void insert (const PositionsVector& positionsVector, long i);
    void append (const PositionsVector& positionsVector);
    void prepend(const PositionsVector& positionsVector);

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override;

    double norm(long i, long j) const;

    const std::vector<PositionsVector>& positionsVectorCollection() const;
    std::vector<PositionsVector>& positionsVectorCollection();

    long numberOfPositionsEntities() const;

private:
    std::vector<PositionsVector> positionsVectorCollection_;
    long numberOfPositionEntities_;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<PositionsVectorCollection> {
        static Node encode(const PositionsVectorCollection &rhs);
        static bool decode(const Node &node, PositionsVectorCollection &rhs);
    };
    Emitter &operator<<(Emitter &out, const PositionsVectorCollection &p) ;
}

#endif //INPSIGHTS_POSITIONSVECTORCOLLECTION_H
