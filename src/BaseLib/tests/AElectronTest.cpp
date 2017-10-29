//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "Electron.h"

using namespace testing;
using namespace Eigen;

class AElectronTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(AElectronTest, Constructor){
    Electron electron(Vector3d(1,2,3),Spin::SpinType::alpha);
}

TEST_F(AElectronTest, CopyConstructor){
    Electron electron(Vector3d(1,2,3),Spin::SpinType::alpha);
    Electron copyElectron(electron);
}

TEST_F(AElectronTest, MemberAccess) {
    Vector3d position = Vector3d(1, 2, 3);
    Spin::SpinType spin = Spin::SpinType::alpha;
    Electron electron(position, spin);

    ASSERT_EQ(electron.spin(), spin);
    ASSERT_EQ(electron.position(), position);
}
