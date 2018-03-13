//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_ELEMENTTYPECOLLECTION_H
#define AMOLQCPP_ELEMENTTYPECOLLECTION_H

#include <Eigen/Core>
#include "AbstractVector.h"
#include "ElementType.h"

class ElementTypesVector : public AbstractVector{
public:
    explicit ElementTypesVector(long size = 0);
    explicit ElementTypesVector(const Eigen::VectorXi& elementTypes);

    Elements::ElementType operator[](long i) const;

    void insert(Elements::ElementType elementType, long i);
    void append(Elements::ElementType elementType);
    void prepend(Elements::ElementType elementType);
    void permute(long i, long j) override ;

    const Eigen::VectorXi& elementTypesAsEigenVector()const;
    Eigen::VectorXi& elementTypesAsEigenVector();


    friend std::ostream& operator<<(std::ostream& os, const ElementTypesVector& ec);

private:
    Eigen::VectorXi elementTypes_;
};

#endif //AMOLQCPP_ELEMENTTYPECOLLECTION_H
