//
// Created by Michael Heuer on 26.09.18.
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


Eigen::PermutationMatrix<Eigen::Dynamic>
BestMatch::swapPermutation(Swap swap, Eigen::Index length) {
    return swapPermutation(swap.first, swap.second, length);
}

Eigen::PermutationMatrix<Eigen::Dynamic>
BestMatch::swapPermutation(Eigen::Index i, Eigen::Index j, Eigen::Index length) {
    assert(i >= 0 && i < length
           && "i must be greater than or equal to zero and smaller than the permutation vector length");
    assert(j >= 0 && j < length
           && "i must be greater than or equal to zero and smaller than the permutation vector length");

    Eigen::PermutationMatrix<Eigen::Dynamic> swapPerm(length);
    swapPerm.setIdentity();
    auto temp = swapPerm.indices()[i];
    swapPerm.indices()[i] = swapPerm.indices()[j];
    swapPerm.indices()[j] = temp;
    return swapPerm;
}

Eigen::PermutationMatrix<Eigen::Dynamic>
BestMatch::concatenateSwaps(std::deque<Swap> swaps, unsigned permutationSize) {

    auto combinedPermutationInKitSystem = Eigen::PermutationMatrix<Eigen::Dynamic>(permutationSize);
    combinedPermutationInKitSystem.setIdentity();

    while (!swaps.empty()) {

        std::pair<Eigen::Index, Eigen::Index> currentSwap = swaps.front();

        assert(currentSwap.first < permutationSize
        && "The first index must be smaller than the total permutation size.");
        assert(currentSwap.second < permutationSize
        && "The second index must be smaller than the total permutation size.");

        // apply current swap
        combinedPermutationInKitSystem =
                combinedPermutationInKitSystem * BestMatch::swapPermutation(currentSwap, permutationSize);

        // pop the current swap from the deque
        swaps.pop_front();

        // treat remaining swaps
        for (auto &s : swaps) {
            if (s.first == currentSwap.first)
                s.first = currentSwap.second;
            else if (s.first == currentSwap.second)
                s.first = currentSwap.first;

            if (s.second == currentSwap.first)
                s.second = currentSwap.second;
            else if (s.second == currentSwap.second)
                s.second = currentSwap.first;

            /* The pairs p_i of a combination have the following properties:
             * p_i = (a,b), a<b
             * p_i+1 = (c,d), a<=b, c<=d
             *
             * Swaps may destroy this order, so we don't make use of it.
             * */
        }
    }
    return combinedPermutationInKitSystem;
}