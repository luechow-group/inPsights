//
// Created by Michael Heuer on 29.05.18.
//

#include <gtest/gtest.h>
#include "MolecularSpectrum.h"

class AMolecularSpectrumTest : public ::testing::Test {
public:

};


TEST_F(AMolecularSpectrumTest , TestName) {

    Particle<Elements::ElementType> H0(Elements::ElementType::H,{1,2,3});
    Particle<Elements::ElementType> H1(Elements::ElementType::H,{4,5,6});

    ExpansionSettings::defaults();

    TypeSpecificNeighborhoodsAtOneCenter tmn;
    tmn[int(Elements::ElementType::H)] = NeighborhoodExpansion();
    tmn[int(Elements::ElementType::He)] = NeighborhoodExpansion();



    NumberedType<int> nt(int(Elements::ElementType::H),0);
    MolecularCenters mn;
    mn[nt] = tmn;

    auto res = mn.find(nt);

    std::cout << res->first << res->second[0] << std::endl;
}
