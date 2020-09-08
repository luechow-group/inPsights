// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <MapUtils.h>

using namespace testing;

TEST(AMapUtilsTest, FindKeysToValueInMap) {
    std::map<int, std::string> map = {
            {0, "b"},
            {1, "a"},
            {2, "b"}};

    auto keysMappedToValueA = MapUtils::findByValue(map, std::string("a"));
    auto keysMappedToValueB = MapUtils::findByValue(map, std::string("b"));
    auto keysMappedToValueC = MapUtils::findByValue(map, std::string("c"));

    ASSERT_THAT(keysMappedToValueA,ElementsAre(1));
    ASSERT_THAT(keysMappedToValueB,ElementsAre(0,2));
    ASSERT_THAT(keysMappedToValueC,IsEmpty());
}
