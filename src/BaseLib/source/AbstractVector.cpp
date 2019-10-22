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

#include <AbstractVector.h>
#include <assert.h>

AbstractVector::AbstractVector(long numberOfEntities, long entityLength)
        :
        entityLength_(entityLength),
        numberOfEntities_(numberOfEntities){
    assert(numberOfEntities >= 0 && "The number of Entities must be non-negative.");
};

void AbstractVector::incrementNumberOfEntities(){
    numberOfEntities_++;
}

long AbstractVector::numberOfEntities() const{
    return numberOfEntities_;
}

long AbstractVector::entityLength() const{
    return entityLength_;
}

void AbstractVector::setEntityLength(long entityLength) {
    assert(entityLength > 0  && "The entity length must be positive.");
    entityLength_ = entityLength;
}

void AbstractVector::setNumberOfEntities(long numberOfEntities){
    assert(numberOfEntities >= 0  && "The number of entities must be non-negative.");
    numberOfEntities_ = numberOfEntities;
}

long AbstractVector::calculateIndex(long i) const {
    assert(i <= numberOfEntities() && "Index is out of bounds");// less or equal because of insert
    assert(i >= -numberOfEntities() && "Reverse index is out of bounds");
    if (i >= 0) return i*entityLength_;
    return (numberOfEntities()+i)*entityLength_;
}

Eigen::PermutationMatrix<Eigen::Dynamic> AbstractVector::randomPermutation(std::default_random_engine& rng) const {
    Eigen::PermutationMatrix<Eigen::Dynamic> perm(numberOfEntities());
    perm.setIdentity();
    std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(),rng);

    return perm;
}
