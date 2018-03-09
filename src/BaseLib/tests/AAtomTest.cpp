//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "Atom.h"
#include "ElementInfo.h"


using namespace testing;
using namespace Eigen;

class AAtomTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(AAtomTest, Charge){
    Atom atom(Vector3d(1,2,3),Elements::ElementType::Ag);
    ASSERT_EQ(Elements::ElementInfo::Z(atom.elementType()),atom.charge());
}
