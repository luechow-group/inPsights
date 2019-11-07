/* Copyright (C) 2017-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
