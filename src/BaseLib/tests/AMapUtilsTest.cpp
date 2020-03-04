/* Copyright (C) 2020 Michael Heuer.
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
