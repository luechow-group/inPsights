/* Copyright (C) 2017-2019 Michael Heuer.
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

#include <gmock/gmock.h>
#include "SpinType.h"

using namespace testing;
using namespace Spins;

class ASpinTypeTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(ASpinTypeTest , IntegerValues){
    ASSERT_EQ(int(SpinType::alpha),-2);
    ASSERT_EQ(int(SpinType::beta),-1);
    ASSERT_EQ(int(SpinType::none),0);
}

TEST_F(ASpinTypeTest , QuantumNumbers){
    ASSERT_EQ(quantumNumber(),1/2.0);
    ASSERT_EQ(magneticQuantumNumber(SpinType::alpha),1/2.0);
    ASSERT_EQ(magneticQuantumNumber(SpinType::beta),-1/2.0);
    EXPECT_DEATH(magneticQuantumNumber(SpinType::none),"");
}
