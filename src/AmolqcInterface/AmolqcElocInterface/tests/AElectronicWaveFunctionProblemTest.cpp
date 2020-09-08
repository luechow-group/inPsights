// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
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

