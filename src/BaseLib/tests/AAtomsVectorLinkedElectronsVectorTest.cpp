//
// Created by Michael Heuer on 26.11.18.
//

#include <gmock/gmock.h>
#include <AtomsVectorLinkedElectronsVector.h>
#include <memory>

using namespace testing;

class AAtomsVectorLinkedElectronsVector: public Test {
public:
    Electron e0,e1,e2;
    ElectronsVector electrons;
    std::shared_ptr<AtomsVector> atomsPtr;

    void SetUp() override {
        electrons = ElectronsVector({
            {Spin::alpha,{1, 0, 0}},
            {Spin::alpha,{1, 0, 0}},
            {Spin::alpha,{1.1, 0, 0}},
            {Spin::beta,{-1, 0, 0}},
            {Spin::beta,{-1.09, 0, 0}}
        });

        atomsPtr = std::make_shared<AtomsVector>(
                AtomsVector({
                    {Element::H ,{1, 0, 0}},
                    {Element::He,{-1, 0, 0}}
                }));

    };
};

TEST_F(AAtomsVectorLinkedElectronsVector, CoreElectrons) {
    double thresh = 0.1;
    auto linkedElectronsVector = AtomsVectorLinkedElectronsVector(atomsPtr,electrons);

    ASSERT_THAT(linkedElectronsVector.coreElectronsIndices(0,thresh), ElementsAre(0,1));
    ASSERT_THAT(linkedElectronsVector.coreElectronsIndices(1,thresh), ElementsAre(3,4));
    ASSERT_THAT(linkedElectronsVector.coreElectronsIndices(thresh), ElementsAre(0,1,3,4));
}

TEST_F(AAtomsVectorLinkedElectronsVector, ValenceElectrons) {
    auto linkedElectronsVector = AtomsVectorLinkedElectronsVector(atomsPtr,electrons);
    ASSERT_THAT(linkedElectronsVector.valenceElectronsIndices(),ElementsAre(2));
}
