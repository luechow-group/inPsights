// Edited by Leonard Reuter on 26.06.19.
// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BestMatch.h>
#include <SOAPSettings.h>

template <>
bool BestMatch::AscendingMetricResult::operator<(const BestMatch::Result<true> &rhs) {
    if(metric == rhs.metric)
        for (Eigen::Index i = 0; i < permutation.indices().size(); ++i) {
            if(permutation.indices()[i] != rhs.permutation.indices()[i])
                return permutation.indices()[i] < rhs.permutation.indices()[i];
        }

    return metric < rhs.metric;
}

template <>
bool BestMatch::DescendingMetricResult::operator<(const BestMatch::Result<false> &rhs) {
    auto numericalPrecisionEpsilon = SOAP::General::settings.comparisonEpsilon.get();

    if(abs(metric - rhs.metric) < numericalPrecisionEpsilon)
        for (Eigen::Index i = 0; i < permutation.indices().size(); ++i) {
            if(permutation.indices()[i] != rhs.permutation.indices()[i])
                return permutation.indices()[i] < rhs.permutation.indices()[i];
        }

    return metric > rhs.metric;
}

Eigen::PermutationMatrix<Eigen::Dynamic> Permutations::combinePermutations(
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

Eigen::PermutationMatrix<Eigen::Dynamic> Permutations::getPermutationToFront(const std::list<long> &relevantIndices, size_t size){
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

Eigen::PermutationMatrix<Eigen::Dynamic> Permutations::getPermutationToBack(const std::list<long> &relevantIndices, size_t size){
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

Eigen::PermutationMatrix<Eigen::Dynamic> Permutations::headToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, size_t size){
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

Eigen::PermutationMatrix<Eigen::Dynamic> Permutations::tailToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, size_t size){
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
