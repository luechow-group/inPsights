// Edited by Leonard Reuter on 26.06.19.
// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_BESTMATCH_H
#define INPSIGHTS_BESTMATCH_H

#include <Eigen/Core>
#include <list>
#include <deque>
#include <ParticlesVector.h>

namespace BestMatch {
    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations( // TODO rename to "concatenatePermutations" ?
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ = false);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    getPermutationToFront(const std::list<long> &relevantIndices, size_t size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    getPermutationToBack(const std::list<long> &relevantIndices, size_t size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    headToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, size_t size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    tailToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, size_t size);

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

    template<bool ascending = true>
    struct Result {
        double metric;
        Eigen::PermutationMatrix<Eigen::Dynamic> permutation;

        bool operator<(const Result &rhs);
    };

    using AscendingMetricResult = Result<true>;
    using DescendingMetricResult = Result<false>;
};

#endif //INPSIGHTS_BESTMATCH_H
