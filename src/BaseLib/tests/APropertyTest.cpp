// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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