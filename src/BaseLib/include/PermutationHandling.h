// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PERMUTATIONHANDLING_H
#define INPSIGHTS_PERMUTATIONHANDLING_H

#include "ParticlesVector.h"
#include <Eigen/Core>
#include <list>

namespace PermutationHandling {
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
}

#endif //INPSIGHTS_PERMUTATIONHANDLING_H
