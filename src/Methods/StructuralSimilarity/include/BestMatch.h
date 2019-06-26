//
// Created by Michael Heuer on 06.09.18.
// Edited by Leonard Reuter on 26.06.19.
//

#ifndef INPSIGHTS_BESTMATCH_H
#define INPSIGHTS_BESTMATCH_H

#include <Eigen/Core>
#include <list>

namespace BestMatch {
    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations(
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ = false);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    getPermutationToFront(const std::list<long> &relevantIndices, const long &size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    getFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const long &size);

    struct Result {
        const double metric;
        const Eigen::PermutationMatrix<Eigen::Dynamic> permutation;
    };


};

#endif //INPSIGHTS_BESTMATCH_H
