/* Copyright (C) 2018-2019 Michael Heuer.
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
#include <Property.h>

TEST(APropertyTest, YamlConversion) {
    Property<int> integer("name");
    integer = 5;

    ASSERT_EQ(integer.get(), 5);
    ASSERT_STREQ(integer.name().c_str(), "name");

    auto node = YAML::convert<Property<int>>::encode(integer);
    Property<int> decodedInteger("name");
    YAML::convert<Property<int>>::decode(node, decodedInteger);

    ASSERT_STREQ(decodedInteger.name().c_str(), integer.name().c_str());
    ASSERT_EQ(decodedInteger.get(), integer.get());
}

enum TestEnum{
    a=0, b=1
};

TEST(APropertyTest, Enum) {
    Property<TestEnum> enumProperty = {TestEnum::a};
    ASSERT_EQ(enumProperty.get(),TestEnum::a);

    enumProperty = TestEnum::b;
    ASSERT_EQ(enumProperty.get(),TestEnum::b);
}