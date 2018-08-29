//
// Created by Michael Heuer on 08.03.17.
//

#ifndef AMOLQCPP_POSITIONSVECTORCOLLECTION_H
#define AMOLQCPP_POSITIONSVECTORCOLLECTION_H

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
    void permute(long i, long j) override;
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

#endif //AMOLQCPP_POSITIONSVECTORCOLLECTION_H
