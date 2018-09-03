//
// Created by Michael Heuer on 08.05.18.
//

#include "TypesVector.h"

template<> SpinTypesVector::TypesVector(
        unsigned long numberOfAlphaElectrons,
        unsigned long numberOfBetaElectrons)
        : IInsertable<int>(0)
{
    for (unsigned long i = 0; i < numberOfAlphaElectrons; ++i)
        this->append(Spin::alpha);

    for (unsigned long i = 0; i < numberOfBetaElectrons; ++i)
        this->append(Spin::beta);
}

template<> SpinTypesVector::TypesVector(std::vector<Spin> types)
        : IInsertable<int>(0)
{
    for (const auto& type : types){
        assert(int(type) >= int(Spins::first()));
        assert(int(type) <= int(Spins::last()));
        this->append(type);
    }
}

template<> ElementTypesVector::TypesVector(std::vector<Element> types)
        : IInsertable<int>(0)
{
    for (const auto& type : types){
        assert(int(type) >= int(Elements::first()));
        assert(int(type) <= int(Elements::last()));
        this->append(type);
    }
}
