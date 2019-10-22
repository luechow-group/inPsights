/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
