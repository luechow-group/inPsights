//
// Created by Michael Heuer on 22.06.18.
//

#include <gtest/gtest.h>
#include <sstream>
#include <TestMolecules.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <sstream>

using namespace testing;

class AYamlTest : public Test {
public:
};

TEST_F(AYamlTest, JsonCompatibility){

    MolecularGeometry mol = TestMolecules::HeH::ElectronsInCores::normal;

    YAML::Emitter out;
    out << YAML::Flow;
    out << YAML::DoubleQuoted;
    out << YAML::BeginMap;
    out << YAML::Key << "Molecule" << YAML::Value << mol;
    out << YAML::EndMap;

    std::string outstring = out.c_str();
    std::string expected = "{\"Molecule\": {\"Atoms\": {\"Types\": [\"He\", \"H\"],\n  \"Positions\": [[0, 0, 0.37], [0, 0, -0.37]]},\n\"Electrons\": {\"Types\": [\"a\", \"a\", \"b\"],\n  \"Positions\": [[0, 0, 0.37], [0, 0, 0.37], [0, 0, -0.37]]}}}";
    ASSERT_EQ(outstring,expected);
}