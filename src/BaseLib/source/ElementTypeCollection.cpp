//
// Created by Michael Heuer on 29.10.17.
//

#include <ElementInfo.h>
#include "Particle.h"
#include "ElementTypeCollection.h"

using namespace Eigen;

ElementTypeCollection::ElementTypeCollection(long size)
        :AbstractCollection(size),
         elementTypes_(VectorXi::Constant(size,int(Elements::ElementType::none)))
{}

ElementTypeCollection::ElementTypeCollection(const VectorXi& elementTypes)
        : AbstractCollection(elementTypes.size()),
          elementTypes_(numberOfEntities_)
{
    assert(elementTypes.minCoeff() >= int(Elements::ElementType::none));
    assert(elementTypes.maxCoeff() <= int(Elements::ElementType::Cn));

    elementTypes_ = elementTypes;
}

Elements::ElementType ElementTypeCollection::elementType(long i) const {
    return  Elements::ElementType(elementTypes_[i]);
}

unsigned long ElementTypeCollection::numberOfElementTypes() const {
    return numberOfEntities_;
}

void ElementTypeCollection::insert(Elements::ElementType elementType, long i) {
    VectorXi before = elementTypes_.head(i);
    VectorXi after = elementTypes_.tail(numberOfEntities_-i);

    elementTypes_.resize(numberOfEntities_+1);
    //elementTypes_ << before, int(elementType), after;
    elementTypes_.head(i) = before;
    elementTypes_.segment(i,1) = Eigen::Matrix<int,1,1>(int(elementType));
    elementTypes_.tail(numberOfEntities_-i) = after;
    ++numberOfEntities_;
}

void ElementTypeCollection::prepend(Elements::ElementType elementType) {
    this->insert(elementType,0);
}

void ElementTypeCollection::append(Elements::ElementType elementType) {
    this->insert(elementType,numberOfEntities_);
}

void ElementTypeCollection::setElementType(long i, Elements::ElementType ElementType) {
    elementTypes_[i] = int(ElementType);
}

VectorXi ElementTypeCollection::elementTypesAsEigenVector() {
    return elementTypes_;
}

void ElementTypeCollection::permute(long i, long j) {
    assert( i >= 0 && i < numberOfEntities_
            && "Index i must be greater than zero and smaller than the number of elements." );
    assert( j >= 0 && j < numberOfEntities_
            && "Index j must be greater than zero and smaller than the number of elements." );
    if(i != j) {
        int temp = elementTypes_[i];
        elementTypes_[i] = elementTypes_[j];
        elementTypes_[j] = temp;
    }
}

std::ostream& operator<<(std::ostream& os, const ElementTypeCollection& ec){
    for (unsigned long i = 0; i < ec.numberOfEntities_; i++) {
        std::string elementSymbol = Elements::ElementInfo::symbol(ec.elementType(i));

        os << elementSymbol
           << std::string(ParticleFormat::significantDigits+3-elementSymbol.length(), ' ')
           << ParticleFormat::separator;
    }
    std::cout << std::endl;
    return os;
}
