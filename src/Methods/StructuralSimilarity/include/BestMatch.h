//
// Created by Michael Heuer on 06.09.18.
//

#ifndef INPSIGHTS_BESTMATCH_H
#define INPSIGHTS_BESTMATCH_H

#include <Eigen/Core>

namespace BestMatch {
    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations(
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ = false);

    struct Result {
        const double metric;
        const Eigen::PermutationMatrix<Eigen::Dynamic> permutation;
    };


};

#endif //INPSIGHTS_BESTMATCH_H
