//
// Created by Leonard Reuter on 12.03.18.
//

#include <gtest/gtest.h>
#include <Eigen/Core>
#include "ElectronicWaveFunctionProblem.h"

using namespace testing;
using namespace Eigen;

TEST(AElectronicWaveFunctionProblemTest, DefaultConstruction) {
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem;

    std::stringstream oss;
    oss << electronicWaveFunctionProblem.getAtomsVector();
    ASSERT_EQ(oss.str(),"");
}

