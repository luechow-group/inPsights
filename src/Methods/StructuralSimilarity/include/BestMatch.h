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
    getPermutationToFront(const std::list<long> &relevantIndices, const long &size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    getPermutationToBack(const std::list<long> &relevantIndices, const long &size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    headToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const long &size);

    Eigen::PermutationMatrix<Eigen::Dynamic>
    tailToFullPermutation(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const long &size);

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

    template <bool ascending=true>
    struct Result {
        double metric;
        Eigen::PermutationMatrix<Eigen::Dynamic> permutation;

        bool operator<(const Result& rhs);
    };

    using AscendingMetricResult = Result<true>;
    using DescendingMetricResult = Result<false>;
};

#endif //INPSIGHTS_BESTMATCH_H
