//
// Created by Michael Heuer on 26.09.18.
// Edited by Leonard Reuter on 26.06.19.
//

#include <BestMatch.h>

bool BestMatch::Result::operator<(const BestMatch::Result &rhs) {
    return metric < rhs.metric;
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
