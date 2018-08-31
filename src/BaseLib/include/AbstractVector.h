//
// Created by Michael Heuer on 05.02.18.
//

#ifndef AMOLQCPP_ABSTRACTVECTOR_H
#define AMOLQCPP_ABSTRACTVECTOR_H

#include <iostream>
#include <Eigen/Core>

/* AbstractVector
 * keeps track of the number of countable entities
 */
class AbstractVector {
public:
    long numberOfEntities() const;

protected:
    explicit AbstractVector(long numberOfEntities = 0);

    void incrementNumberOfEntities();

    void setNumberOfEntities(long numberOfEntities);

    virtual void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& permutation) = 0;

    virtual long calculateIndex(long i) const;

private:
    long numberOfEntities_;
};

#endif //AMOLQCPP_ABSTRACTVECTOR_H
