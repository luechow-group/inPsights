//
// Created by Michael Heuer on 2018-12-22.
//

#include <gmock/gmock.h>
#include <GeneralSettings.h>
#include <Property.h>

class TestSettings : public GeneralSettings {
public:
    Property<int> number = {1234567890, VARNAME(number)};
    Property<double> threshold = {1.234567890, VARNAME(threshold)};

    TestSettings() = default;

    explicit TestSettings(const YAML::Node &node) {
        intProperty::decode(node, number);
        doubleProperty::decode(node, threshold);
    }

    void addToNode(YAML::Node &node) const override{
        node[number.name()] = number.get();
        node[threshold.name()] = threshold.get();
    }
};
YAML_GENERALSETTINGS_DECLARATION(TestSettings)
YAML_GENERALSETTINGS_DEFINITION(TestSettings)


TEST(AGeneralSettingsTest, YamlConversion) {
    TestSettings settings;
    ASSERT_STREQ(settings.number.name().c_str(), "number");
    ASSERT_EQ(settings.number.get(), 1234567890);
    ASSERT_STREQ(settings.threshold.name().c_str(), "threshold");
    ASSERT_EQ(settings.threshold.get(),1.234567890);

    settings.number = 123;
    settings.threshold = -1.23;
    auto node = YAML::convert<TestSettings>::encode(settings);

    ASSERT_TRUE(node["number"]);
    ASSERT_EQ(node["number"].as<int>(), settings.number.get());
    ASSERT_TRUE(node["threshold"]);
    ASSERT_EQ(node["threshold"].as<double>(), settings.threshold.get());

    TestSettings decodedSettings(node);

    ASSERT_STREQ(decodedSettings.number.name().c_str(), settings.number.name().c_str());
    ASSERT_EQ(decodedSettings.number.get(), settings.number.get());
    ASSERT_STREQ(decodedSettings.threshold.name().c_str(), settings.threshold.name().c_str());
    ASSERT_EQ(decodedSettings.threshold.get(), settings.threshold.get());
}
