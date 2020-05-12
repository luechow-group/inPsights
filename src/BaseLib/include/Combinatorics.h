/* Copyright (C) 2019 Michael Heuer.
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

#ifndef INPSIGHTS_COMBINATORICS_H
#define INPSIGHTS_COMBINATORICS_H

#include <vector>
#include <Eigen/Core>
#include <algorithm>

namespace Combinatorics {
    const size_t MAX_BINOMIAL = 34;
    const size_t MAX_FACTORIAL = 12;

    // C(n, k) = n! / (n-k)! / !k
    std::size_t binomial(std::size_t n, std::size_t k);
    std::size_t factorial(std::size_t n);

    template<
            typename ObjectType,
            typename SubcontainerType = std::vector<ObjectType>,
            typename ContainerType = std::vector<SubcontainerType>> class Combinations {
    public:
        Combinations(const SubcontainerType& objects, unsigned numberOfElements)
        :
        N(objects.size()),
        M(numberOfElements),
        count(0),
        set(objects),
        partial(SubcontainerType(M)),
        out(binomial(N, M), SubcontainerType(M)) {
            assert(M <= N && "The combinations cannot include more elements than the maximal number of elements.");

            if(M > 0)
                generate(0, N-1, M-1);
            //else: the out vector contains one element with a size of zero
        };

        using const_iterator = typename ContainerType::const_iterator ;
        using iterator = typename ContainerType::iterator;
        iterator begin() { return out.begin(); }
        iterator end() { return out.end(); }
        const ContainerType & get() const { return out; }
        ContainerType & get() { return out; }

    private:
        void generate(std::size_t i, std::size_t j, std::size_t m) {
            // combination of size m (number of slots) out of set[i..j]
            if (m > 0) {
                for (std::size_t z=i; z<j-m+1; z++) {
                    partial[M-m-1]=set[z]; // add element to permutation
                    generate(z+1, j, m-1);
                }
            } else {
                // last position
                for (std::size_t z=i; z<j-m+1; z++) {
                    partial[M-m-1] = set[z];
                    out[count++] = SubcontainerType(partial); // add to output vector
                }
            }
        }

        std::size_t N,M,count;
        SubcontainerType set,partial;
        ContainerType out;
    };

    template<
            typename ObjectType,
            typename SubcontainerType = std::vector<ObjectType>,
            typename ContainerType = std::vector<SubcontainerType>> class Permutations {
    public:
        explicit Permutations(const SubcontainerType& objects):
                N(objects.size()),
                M(factorial(N)),
                set(objects),
                out(factorial(N), SubcontainerType(N)) {
            generate();
        };

        using const_iterator = typename ContainerType::const_iterator ;
        using iterator = typename ContainerType::iterator;
        iterator begin() { return out.begin(); }
        iterator end() { return out.end(); }
        const ContainerType & get() const { return out; }
        ContainerType & get() { return out; }

    private:
        void generate(){
            for (std::size_t i = 0; i < M; ++i) {
                out[i] = set;
                std::next_permutation(set.begin(), set.end());
            }
        };

        std::size_t N,M;
        SubcontainerType set;
        ContainerType out;
    };
}

#endif //INPSIGHTS_COMBINATORICS_H
