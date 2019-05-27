//
// Created by Michael Heuer on 06.09.18.
//

#ifndef INPSIGHTS_BESTMATCH_H
#define INPSIGHTS_BESTMATCH_H

#include <Eigen/Core>
#include <deque>
#include <ParticlesVector.h>

namespace BestMatch {
    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations( // TODO rename to "concatenatePermutations" ?
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ = false);

    using Swap = std::pair<Eigen::Index, Eigen::Index>;

    Eigen::PermutationMatrix<Eigen::Dynamic> swapPermutation(Swap swap, Eigen::Index length);
    Eigen::PermutationMatrix<Eigen::Dynamic> swapPermutation(Eigen::Index i, Eigen::Index j, Eigen::Index length);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    concatenateSwaps(std::deque<Swap> swaps, unsigned permutationSize);


    template<typename Type>
    Eigen::PermutationMatrix<Eigen::Dynamic> findTypeSeparatingPermutation(
            const ParticlesVector<Type> &particlesVector) {

        Eigen::VectorXi typeSerparatingPermutationIndices(particlesVector.numberOfEntities());

        Eigen::Index i = 0;

        for (const auto&[type, count] : particlesVector.typesVector().countTypes()) {

            for (std::size_t j = 0; j < count; ++j) {
                auto[foundQ, index] = particlesVector.typesVector().findIndexOfEnumeratedType(
                        EnumeratedType<Type>(type, j));
                assert(foundQ);
                typeSerparatingPermutationIndices[i] = index;
                ++i;
            }
        }
        return Eigen::PermutationMatrix<Eigen::Dynamic>(typeSerparatingPermutationIndices);
    }

    struct Result {
        double metric;
        Eigen::PermutationMatrix<Eigen::Dynamic> permutation;

        bool operator<(const Result& rhs);

    };


};

#endif //INPSIGHTS_BESTMATCH_H
