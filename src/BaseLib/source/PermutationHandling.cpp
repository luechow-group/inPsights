// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PermutationHandling.h"

Eigen::PermutationMatrix<Eigen::Dynamic> PermutationHandling::combinePermutations(
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

Eigen::PermutationMatrix<Eigen::Dynamic> PermutationHandling::getPermutationToFront(const std::list<long> &relevantIndices, size_t size){
    assert(relevantIndices.size() <= size);

    Eigen::VectorXi indices(size);
    size_t i = 0;
    for (auto index : relevantIndices){
        indices[i] = index;
        i++;
    }

    for (size_t j=0; j<size; j++) {
        if (std::find(relevantIndices.begin(), relevantIndices.end(), j) == relevantIndices.end()){
            indices[i] = j;
            i++;
        }
    }

    assert(i == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices).inverse();
};

Eigen::PermutationMatrix<Eigen::Dynamic> PermutationHandling::getPermutationToBack(const std::list<long> &relevantIndices, size_t size){
    assert(relevantIndices.size() <= size);

    Eigen::VectorXi indices(size);
    size_t i = 0;

    for (size_t j=0; j<size; j++) {
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

Eigen::PermutationMatrix<Eigen::Dynamic> PermutationHandling::headToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, size_t size){
    // builds full permutation from head permutation
    assert(size_t(permutation.indices().size()) <= size);

    Eigen::VectorXi indices(size);
    long i = 0;
    while (i<permutation.indices().size()){
        indices[i] = permutation.indices()[i];
        i++;
    }

    for (size_t j=permutation.indices().size(); j<size; j++) {
        indices[i] = j;
        i++;
    }

    assert(size_t(i) == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
};

Eigen::PermutationMatrix<Eigen::Dynamic> PermutationHandling::tailToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, size_t size){
    // builds full permutation from tail permutation
    assert(size_t(permutation.indices().size()) <= size);

    Eigen::VectorXi indices(size);

    long headSize = size - permutation.indices().size();
    long i = 0;
    for (long j=0; j < headSize; j++) {
        indices[i] = i;
        i++;
    }

    for (long j=0; j < permutation.indices().size(); j++) {
        indices[i] = permutation.indices()[j] + headSize;
        i++;
    }

    assert(size_t(i) == size);

    return Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
};
