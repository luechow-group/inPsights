/* Copyright (C) 2019 Michael Heuer.
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
#include <Enumerate.h>
#include <vector>

struct Thing{
    double number;
    std::string name;

    friend bool operator==(const Thing& lhs, const Thing& rhs) {
        return std::tie(lhs.number, lhs.name) == std::tie(rhs.number, rhs.name);
    }
};

TEST(AEnumerateTest, Enumerate) {
    std::vector<Thing> things =  {{1.2, "Joe"},{-2.6, "Bart"}};

    for(auto [i, thing] : enumerate(things))
        ASSERT_EQ(thing, things[i]);
}