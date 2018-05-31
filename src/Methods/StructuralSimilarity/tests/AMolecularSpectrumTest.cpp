//
// Created by Michael Heuer on 29.05.18.
//

#include <gtest/gtest.h>
#include "MolecularSpectrum.h"

class AMolecularSpectrumTest : public ::testing::Test {
public:

};


TEST_F(AMolecularSpectrumTest , TestName) {

    Atom H0(Element::H,{1,2,3});
    Atom H1(Element::H,{4,5,6});

    ExpansionSettings::defaults();

    TypeSpecificNeighborhoodsAtOneCenter tmn;
    tmn[int(Element::H)] = NeighborhoodExpansion();
    tmn[int(Element::He)] = NeighborhoodExpansion();



    NumberedType<int> nt(int(Element::H),0);
    MolecularCenters mn;
    mn[nt] = tmn;

    auto res = mn.find(nt);

    std::cout << res->first << res->second[0] << std::endl;
}
