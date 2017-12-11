//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELEMENTTYPECOLLECTION_H
#define AMOLQCGUI_ELEMENTTYPECOLLECTION_H

#include <Eigen/Core>
#include "ElementType.h"

class ElementTypeCollection{
public:
    explicit ElementTypeCollection(long size = 0);
    explicit ElementTypeCollection(const Eigen::VectorXi& elementTypes);

    Elements::ElementType elementType(long i) const;

    void insert(Elements::ElementType elementType, long i);
    void append(Elements::ElementType elementType);
    void prepend(Elements::ElementType elementType);
    void permute(long i, long j);
    
    void setElementType(long i, Elements::ElementType elementType);

    Eigen::VectorXi elementTypesAsEigenVector();

    unsigned long numberOfElementTypes() const;

private:
    unsigned long numberOfElementsTypes_;
    Eigen::VectorXi elementTypes_;
};

#endif //AMOLQCGUI_ELEMENTTYPECOLLECTION_H
