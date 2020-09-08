// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <TypesVector.h>

template<> SpinTypesVector::TypesVector(
        unsigned long numberOfAlphaElectrons,
        unsigned long numberOfBetaElectrons)
        : InsertableVector<int>(0)
{
    for (unsigned long i = 0; i < numberOfAlphaElectrons; ++i)
        this->append(Spin::alpha);

    for (unsigned long i = 0; i < numberOfBetaElectrons; ++i)
        this->append(Spin::beta);
}

template<> SpinTypesVector::TypesVector(std::vector<Spin> types)
        : InsertableVector<int>(0)
{
    for (const auto& type : types){
        assert(int(type) >= int(Spins::first()));
        assert(int(type) <= int(Spins::last()));
        this->append(type);
    }
}

template<>
unsigned SpinTypesVector::multiplicity() {
    auto numberOfUnpairedElectrons = unsigned(std::abs(
            int(countOccurence(Spin::alpha)) - int(countOccurence(Spin::beta))
            ));
    return numberOfUnpairedElectrons + 1;
}

template<>
void SpinTypesVector::flipSpins() {
    for (long i = 0; i < numberOfEntities(); ++i)
        switch(type(i)) {
            case Spin::alpha: {
                data_[i] = Spins::spinToInt(Spin::beta);
                break;
            }
            case Spin::beta: {
                data_[i] = Spins::spinToInt(Spin::alpha);
                break;
            }
            default: break;
        }
}

template<> ElementTypesVector::TypesVector(std::vector<Element> types)
        : InsertableVector<int>(0)
{
    for (const auto& type : types){
        assert(int(type) >= int(Elements::first()));
        assert(int(type) <= int(Elements::last()));
        this->append(type);
    }
}
