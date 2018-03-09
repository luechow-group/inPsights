//
// Created by Michael Heuer on 29.10.17.
//

#include "ElementTypeCollection.h"
#include "ElementInfo.h"

using namespace Eigen;

ElementTypeCollection::ElementTypeCollection(long numberOfEntities)
        : AbstractCollection(numberOfEntities),
          elementTypes_(VectorXi::Constant(numberOfEntities, int(Elements::ElementType::none)))
{}

ElementTypeCollection::ElementTypeCollection(const VectorXi& elementTypes)
        : AbstractCollection(elementTypes.size()),
          elementTypes_(elementTypes)
{
    assert(elementTypes_.minCoeff() >= int(Elements::first()));
    assert(elementTypes_.maxCoeff() <= int(Elements::last()));
}

Elements::ElementType ElementTypeCollection::operator[](long i) const {
    return  Elements::ElementType(elementTypes_[calculateIndex(i)]);
}

void ElementTypeCollection::insert(Elements::ElementType elementType, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    VectorXi before = elementTypes_.head(i);
    VectorXi after = elementTypes_.tail(numberOfEntities()-i);

    elementTypes_.resize(numberOfEntities()+1);
    elementTypes_ << before, int(elementType), after;

    incrementNumberOfEntities();
}

void ElementTypeCollection::prepend(Elements::ElementType elementType) {
    this->insert(elementType,0);
}

void ElementTypeCollection::append(Elements::ElementType elementType) {
    this->insert(elementType,numberOfEntities());
}

const VectorXi& ElementTypeCollection::elementTypesAsEigenVector() const {
    return elementTypes_;
}

VectorXi& ElementTypeCollection::elementTypesAsEigenVector() {
    return elementTypes_;
}

void ElementTypeCollection::permute(long i, long j) {
    if(i != j) {
        int temp = elementTypes_[calculateIndex(i)];
        elementTypes_[calculateIndex(i)] = elementTypes_[calculateIndex(j)];
        elementTypes_[calculateIndex(j)] = temp;
    }
}

std::ostream& operator<<(std::ostream& os, const ElementTypeCollection& etc){
    std::string symbol;
    for (unsigned long i = 0; i < etc.numberOfEntities(); i++) {
        if (etc[i] != Elements::ElementType::none){
            symbol = Elements::ElementInfo::symbol(etc[i]);
        }
        os << symbol << std::endl;
    }
    return os;
}
