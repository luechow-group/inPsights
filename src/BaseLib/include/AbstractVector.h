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

#ifndef INPSIGHTS_ABSTRACTVECTOR_H
#define INPSIGHTS_ABSTRACTVECTOR_H

#include <iostream>
#include <Eigen/Core>
#include <random>

/* AbstractVector
 * keeps track of the number of countable entities
 */
class AbstractVector {
public:
    long numberOfEntities() const;
    long entityLength() const;
    Eigen::PermutationMatrix<Eigen::Dynamic> randomPermutation(std::default_random_engine& rng) const;

protected:
    explicit AbstractVector(long numberOfEntities = 0, long entityLength = 1);

    void incrementNumberOfEntities();

    void setEntityLength(long entityLength);
    void setNumberOfEntities(long numberOfEntities);

    virtual void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& permutation) = 0;

    long calculateIndex(long i) const;

private:
    long entityLength_;
    long numberOfEntities_;
};

#endif //INPSIGHTS_ABSTRACTVECTOR_H
