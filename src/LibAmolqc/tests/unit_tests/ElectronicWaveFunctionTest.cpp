//
// Created by Michael Heuer on 09.01.17.
//

#include <gmock/gmock.h>
#include "ElectronicWaveFunction.h"

using namespace testing;

class ElectronicWaveFunctionTest : public Test {
public:
  ElectronicWaveFunction wf;
  void SetUp() override {
    wf = ElectronicWaveFunction();
  }
};

TEST_F(ElectronicWaveFunctionTest, IsInitialized) {
  ASSERT_TRUE(false);
}


TEST_F(ElectronicWaveFunctionTest, Psi) {
  ASSERT_TRUE(false);
}


