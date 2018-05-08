//
// Created by Michael Heuer on 08.05.18.
//

#include "TypesVector.h"

template<> SpinTypesVector::TypesVector(unsigned long numberOfAlphaElectrons,
                             unsigned long numberOfBetaElectrons)
        : AbstractVector(0),
          types_(Eigen::VectorXi::Constant(0, 0))
{
    for (unsigned long i = 0; i < numberOfAlphaElectrons; ++i)
        this->append(Spins::SpinType::alpha);

    for (unsigned long i = 0; i < numberOfBetaElectrons; ++i)
        this->append(Spins::SpinType::beta);
}

template<> SpinTypesVector::TypesVector(std::vector<Spins::SpinType> types)
        : AbstractVector(0),
          types_(0)
{
    for (const auto& type : types){
        assert(int(type) >= int(Spins::first()));
        assert(int(type) <= int(Spins::last()));
        this->append(type);
    }
}

template<> ElementTypesVector ::TypesVector(std::vector<Elements::ElementType> types)
        : AbstractVector(0),
          types_(0)
{
    for (const auto& type : types){
        assert(int(type) >= int(Elements::first()));
        assert(int(type) <= int(Elements::last()));
        this->append(type);
    }
}
