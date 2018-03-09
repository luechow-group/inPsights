//
// Created by Michael Heuer on 05.02.18.
//

#ifndef AMOLQCPP_ABSTRACTCOLLECTION_H
#define AMOLQCPP_ABSTRACTCOLLECTION_H

#include <iostream>

/* AbstractVector
 * keeps track of the number of countable entities
 */
class AbstractVector {
public:
    explicit AbstractVector(long numberOfEntities = 0);

    long numberOfEntities() const;

protected:
    void incrementNumberOfEntities();

    void setNumberOfEntities(long numberOfEntities);

    virtual void permute(long i, long j){}; //TODO MAKE IT ABSTRACT

    virtual long calculateIndex(long i) const;

private:
    long numberOfEntities_;
};

#endif //AMOLQCPP_ABSTRACTCOLLECTION_H
