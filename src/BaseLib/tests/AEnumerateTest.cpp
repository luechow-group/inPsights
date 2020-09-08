// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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