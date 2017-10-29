//
// Created by Michael Heuer on 29.10.17.
//

#include "ElementTypeCollection.h"

using namespace Eigen;

ElementTypeCollection::ElementTypeCollection(long size)
        :size_(size),
         elementTypes_(VectorXi::Constant(size,int(Elements::ElementType::none)))
{}

ElementTypeCollection::ElementTypeCollection(const VectorXi& elementTypes)
        : size_(elementTypes.size()),
          elementTypes_(size_)
{
    assert(elementTypes.minCoeff() >= int(Elements::ElementType::none));
    assert(elementTypes.maxCoeff() <= int(Elements::ElementType::Cn));

    elementTypes_ = elementTypes;
}

Elements::ElementType ElementTypeCollection::elementType(long i) {
    return  Elements::ElementType(elementTypes_[i]);
}

void ElementTypeCollection::setElementType(long i, Elements::ElementType ElementType) {
    elementTypes_[i] = int(ElementType);
}

VectorXi ElementTypeCollection::asVectorXi() {
    return elementTypes_;
}