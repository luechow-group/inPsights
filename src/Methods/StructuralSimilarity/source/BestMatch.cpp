// Edited by Leonard Reuter on 26.06.19.
/* Copyright (C) 2018-2019 Michael Heuer.
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

#include <BestMatch.h>

template <>
bool BestMatch::AscendingMetricResult::operator<(const BestMatch::Result<true> &rhs) {

    if(metric != rhs.metric)
        return metric < rhs.metric;
    else
        for (Eigen::Index i = 0; i < permutation.indices().size(); ++i) {
            if(permutation.indices()[i] != rhs.permutation.indices()[i])
                return permutation.indices()[i] < rhs.permutation.indices()[i];
        }

    return metric < rhs.metric;
}

template <>
bool BestMatch::DescendingMetricResult::operator<(const BestMatch::Result<false> &rhs) {

    if(metric != rhs.metric)
        return metric > rhs.metric;
    else
        for (Eigen::Index i = 0; i < permutation.indices().size(); ++i) {
            if(permutation.indices()[i] != rhs.permutation.indices()[i])
                return permutation.indices()[i] < rhs.permutation.indices()[i];
        }

    return metric > rhs.metric;
}

Eigen::PermutationMatrix<Eigen::Dynamic> BestMatch::combinePermutations(
        const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ) {

    long n1 = p1.size(), n2 = p2.size();
    Eigen::VectorXi combined(n1 + n2);

    if (!flipSpinsQ) {
        combined.segment(0, n1) = p1.indices().base();
        combined.segment(n1, n2) = (p2.indices().base().array() + n1);
    } else {
        combined.segment(0, n1) = (p1.indices().base().array() + n1);
        combined.segment(n1, n2) = p2.indices().base();
    }
    return Eigen::PermutationMatrix<Eigen::Dynamic>(combined);
};

Eigen::PermutationMatrix<Eigen::Dynamic> BestMatch::getPermutationToFront(const std::list<long> &relevantIndices, const long &size){
    assert(relevantIndices.size() <= size);

    Eigen::VectorXi indices(size);
    int i = 0;
    for (auto index : relevantIndices){
        indices[i] = index;
        i++;
    }

    for (long j=0; j<size; j++) {
        if (std::find(relevantIndices.begin(), relevantIndices.end(), j) == relevantIndices.end()){
            indices[i] = j;
            i++;
        }
    }

    assert(i == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices).inverse();
};

Eigen::PermutationMatrix<Eigen::Dynamic> BestMatch::getPermutationToBack(const std::list<long> &relevantIndices, const long &size){
    assert(relevantIndices.size() <= size);

    Eigen::VectorXi indices(size);
    int i = 0;

    for (long j=0; j<size; j++) {
        if (std::find(relevantIndices.begin(), relevantIndices.end(), j) == relevantIndices.end()){
            indices[i] = j;
            i++;
        }
    }

    for (auto index : relevantIndices){
        indices[i] = index;
        i++;
    }

    assert(i == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices).inverse();
};

Eigen::PermutationMatrix<Eigen::Dynamic> BestMatch::headToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const long &size){
    // builds full permutation from head permutation
    assert(permutation.indices().size() <= size);

    Eigen::VectorXi indices(size);
    int i = 0;
    while (i<permutation.indices().size()){
        indices[i] = permutation.indices()[i];
        i++;
    }

    for (long j=permutation.indices().size(); j<size; j++) {
        indices[i] = j;
        i++;
    }

    assert(i == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
};

Eigen::PermutationMatrix<Eigen::Dynamic> BestMatch::tailToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const long &size){
    // builds full permutation from tail permutation
    assert(permutation.indices().size() <= size);

    Eigen::VectorXi indices(size);

    int headSize = size - permutation.indices().size();
    int i = 0;
    for (int j=0; j < headSize; j++) {
        indices[i] = i;
        i++;
    }

    for (int j=0; j < permutation.indices().size(); j++) {
        indices[i] = permutation.indices()[j] + headSize;
        i++;
    }

    assert(i == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
};
