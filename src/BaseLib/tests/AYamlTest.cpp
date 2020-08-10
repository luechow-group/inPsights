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
#include <sstream>
#include <TestMolecules.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <sstream>

using namespace testing;

TEST(AYamlTest, JsonCompatibility){

    MolecularGeometry mol = TestMolecules::HeH::ElectronsInCores::normal;

    YAML::Emitter out;
    out << YAML::Flow;
    out << YAML::DoubleQuoted;
    out << YAML::BeginMap;
    out << YAML::Key << "Molecule" << YAML::Value << mol;
    out << YAML::EndMap;

    std::string outstring = out.c_str();
    std::string expected ="{\"Molecule\": {\"Atoms\": {\"Types\": [\"He\", \"H\"],\n  \"Positions\": [\n      [0, 0, 0.37],\n      [0, 0, -0.37],\n      ]},\n\"Electrons\": {\"Types\": [\"a\", \"a\", \"b\"],\n  \"Positions\": [\n      [0, 0, 0.37],\n      [0, 0, 0.37],\n      [0, 0, -0.37],\n      ]}}}";
    ASSERT_EQ(outstring,expected);
}