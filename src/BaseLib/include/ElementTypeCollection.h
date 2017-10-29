//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELEMENTTYPECOLLECTION_H
#define AMOLQCGUI_ELEMENTTYPECOLLECTION_H

#include <Eigen/Core>
#include "ElementType.h"

using namespace Eigen;

class ElementTypeCollection{
public:
    explicit ElementTypeCollection(long size);

    explicit ElementTypeCollection(const VectorXi& spinTypes);

    Elements::ElementType elementType(long i);

    void setElementType(long i, Elements::ElementType elementType);

    VectorXi asVectorXi();

private:
    long size_;
    VectorXi elementTypes_;
};

#endif //AMOLQCGUI_ELEMENTTYPECOLLECTION_H
