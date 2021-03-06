// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <ISettings.h>
#include <Property.h>

namespace Settings {
    class TestSettings : public ISettings {
    public:
        Property<int> number = {1234567890, VARNAME(number)};
        Property<double> threshold = {1.234567890, VARNAME(threshold)};

        void onThresholdChanged(double val){
            assert(val > 0 && "The threshold cannot be negative.");
        };

        TestSettings()
        : ISettings(VARNAME(TestSettings)) {
            threshold.onChange_.connect(std::bind(&TestSettings::onThresholdChanged, this, std::placeholders::_1));
        };

        explicit TestSettings(const YAML::Node &node)
        : TestSettings() {
            intProperty::decode(node[className], number);
            doubleProperty::decode(node[className], threshold);
        };

        void appendToNode(YAML::Node &node) const override {
            node[className][number.name()] = number.get();
            node[className][threshold.name()] = threshold.get();
        }
    };
}
YAML_SETTINGS_DECLARATION(Settings::TestSettings)
YAML_SETTINGS_DEFINITION(Settings::TestSettings)

class TestMethod{
public:
    static inline Settings::TestSettings settings {};
};

TEST(AISettingsTest, YamlConversion) {
    using namespace Settings;
    TestSettings settings;
    ASSERT_STREQ(settings.number.name().c_str(), "number");
    ASSERT_EQ(settings.number.get(), 1234567890);
    ASSERT_STREQ(settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(settings.threshold.get(),1.234567890);

    settings.number = 123;
    settings.threshold = 1.23;
    auto node = YAML::convert<TestSettings>::encode(settings);


    ASSERT_TRUE(node[TestMethod::settings.name()][TestMethod::settings.number.name()]);
    ASSERT_EQ(node[TestMethod::settings.name()][TestMethod::settings.number.name()].as<int>(), settings.number.get());
    ASSERT_TRUE(node[TestMethod::settings.name()][TestMethod::settings.threshold.name()]);
    ASSERT_EQ(node[TestMethod::settings.name()][TestMethod::settings.threshold.name()].as<double>(), settings.threshold.get());

    TestSettings decodedSettings(node);

    ASSERT_STREQ(decodedSettings.number.name().c_str(), settings.number.name().c_str());
    ASSERT_EQ(decodedSettings.number.get(), settings.number.get());
    ASSERT_STREQ(decodedSettings.threshold.name().c_str(), settings.threshold.name().c_str());
    ASSERT_EQ(decodedSettings.threshold.get(), settings.threshold.get());

    EXPECT_DEATH(settings.threshold = -0.1, "The threshold cannot be negative.");
}

TEST(AISettingsTest, StaticMembership) {
    ASSERT_STREQ(TestMethod::settings.number.name().c_str(), "number");
    ASSERT_EQ(TestMethod::settings.number.get(), 1234567890);
    ASSERT_STREQ(TestMethod::settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(TestMethod::settings.threshold.get(), 1.234567890);

    TestMethod::settings.number = 123;
    TestMethod::settings.threshold = 1.23;

    ASSERT_STREQ(TestMethod::settings.number.name().c_str(), "number");
    ASSERT_EQ(TestMethod::settings.number.get(), 123);
    ASSERT_STREQ(TestMethod::settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(TestMethod::settings.threshold.get(), 1.23);
}
