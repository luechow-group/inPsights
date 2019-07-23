//
// Created by Michael Heuer on 05.02.18.
//

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
