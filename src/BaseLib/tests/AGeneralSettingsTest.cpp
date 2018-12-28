//
// Created by Michael Heuer on 2018-12-22.
//

#include <gmock/gmock.h>
#include <ISettings.h>
#include <Property.h>

namespace Settings {
    class TestSettings : public ISettings {
    public:
        inline static const std::string className = {VARNAME(TestSettings)};
        Property<int> number = {1234567890, VARNAME(number)};
        Property<double> threshold = {1.234567890, VARNAME(threshold)};

        TestSettings() {
            threshold.onChange().connect([&](double val) {
                assert(val > 0 && "The threshold cannot be negative.");
            });
        }

        explicit TestSettings(const YAML::Node &node) {
            intProperty::decode(node[className], number);
            doubleProperty::decode(node[className], threshold);
        }

        void addToNode(YAML::Node &node) const override {
            node[className][number.name()] = number.get();
            node[className][threshold.name()] = threshold.get();
        }
    };
}
YAML_SETTINGS_DECLARATION(Settings::TestSettings)
YAML_SETTINGS_DEFINITION(Settings::TestSettings)

class TestMethod{
public:
    static inline Settings::TestSettings settings;
};

TEST(AGeneralSettingsTest, YamlConversion) {
    using namespace Settings;
    TestSettings settings;
    ASSERT_STREQ(settings.number.name().c_str(), "number");
    ASSERT_EQ(settings.number.get(), 1234567890);
    ASSERT_STREQ(settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(settings.threshold.get(),1.234567890);

    settings.number = 123;
    settings.threshold = 1.23;
    auto node = YAML::convert<TestSettings>::encode(settings);


    ASSERT_TRUE(node[TestSettings::className]["number"]);
    ASSERT_EQ(node[TestSettings::className]["number"].as<int>(), settings.number.get());
    ASSERT_TRUE(node[TestSettings::className]["threshold"]);
    ASSERT_EQ(node[TestSettings::className]["threshold"].as<double>(), settings.threshold.get());

    TestSettings decodedSettings(node);

    ASSERT_STREQ(decodedSettings.number.name().c_str(), settings.number.name().c_str());
    ASSERT_EQ(decodedSettings.number.get(), settings.number.get());
    ASSERT_STREQ(decodedSettings.threshold.name().c_str(), settings.threshold.name().c_str());
    ASSERT_EQ(decodedSettings.threshold.get(), settings.threshold.get());

    EXPECT_DEATH(settings.threshold = -0.1, "The threshold cannot be negative.");
}

TEST(AGeneralSettingsTest, StaticMembership) {
    using namespace Settings;
    TestSettings& settings = TestMethod::settings;

    ASSERT_STREQ(settings.number.name().c_str(), "number");
    ASSERT_EQ(settings.number.get(), 1234567890);
    ASSERT_STREQ(settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(settings.threshold.get(), 1.234567890);

    settings.number = 123;
    settings.threshold = 1.23;

    ASSERT_STREQ(settings.number.name().c_str(), "number");
    ASSERT_EQ(settings.number.get(), 123);
    ASSERT_STREQ(settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(settings.threshold.get(), 1.23);
}
