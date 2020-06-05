/* Copyright (C) 2019-2020 Michael Heuer.
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
#include <spdlog/spdlog.h>
#include <Group.h>
#include <Reference.h>
#include <TestMolecules.h>

using namespace testing;

class AGroupTest : public ::testing::Test {
public:
    Group g1,g2,g3, h1,h2,h3;
    AtomsVector atoms;

    void SetUp() override {
        g1 = Group({atoms, 1.1,TestMolecules::H2::ElectronsInCores::normal.electrons(), 0});
        g2 = Group({atoms, 1.2,TestMolecules::H2::ElectronsInCores::translated.electrons(),1});
        g3 = Group({atoms, 1.0,TestMolecules::H2::ElectronsInCores::flippedSpins.electrons(), 2});
        h1 = Group({g1});
        h2 = Group({g2});
        h3 = Group({g3});
    }
};


TEST_F(AGroupTest, DefaultConstructor) {
    Group supergroup0;
    ASSERT_EQ(supergroup0.representative(), nullptr);
    ASSERT_EQ(supergroup0.getSelectedElectronsCount(), 0);
    ASSERT_TRUE(supergroup0.isLeaf());
}

TEST_F(AGroupTest, CopyConstructor){
    Group group(g1);
    ASSERT_TRUE(group.isLeaf());
    ASSERT_EQ(group.representative()->value(), 1.1);
    ASSERT_EQ(group.getSelectedElectronsCount(), 2);
}

TEST_F(AGroupTest, Preallocation){
    Group supergroup(3);

    ASSERT_FALSE(supergroup.isLeaf());
    ASSERT_TRUE(supergroup[0].isLeaf());
    ASSERT_EQ(supergroup.getSelectedElectronsCount(), 0);

    ASSERT_EQ(supergroup.representative(), supergroup.front().representative());
    ASSERT_EQ(supergroup.front().representative(), nullptr);
    ASSERT_EQ(supergroup[1].representative(), nullptr);
    ASSERT_EQ(supergroup[2].representative(), nullptr);

    supergroup[0] = g1;
    supergroup[1] = g2;
    supergroup[2] = g3;

    ASSERT_FALSE(supergroup.isLeaf());
    ASSERT_EQ(supergroup.representative(), g1.representative());
    ASSERT_EQ(supergroup[0].representative(), g1.representative());
    ASSERT_EQ(supergroup.front().getSelectedElectronsCount(), 2);
}

TEST_F(AGroupTest, ListInitialization){
    Group supergroup = {g1, g2, g3};
    ASSERT_FALSE(supergroup.isLeaf());

    ASSERT_EQ(supergroup.representative()->value(), 1.1);
    ASSERT_EQ(supergroup.getSelectedElectronsCount(), 2);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[1].representative()->value(), 1.2);
    ASSERT_EQ(supergroup[2].representative()->value(), 1.0);
}

TEST_F(AGroupTest, ListInitialization_Constructor){
    Group supergroup({g1, g2, g3});
    ASSERT_FALSE(supergroup.isLeaf());

    ASSERT_EQ(supergroup.representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[1].representative()->value(), 1.2);
    ASSERT_EQ(supergroup[2].representative()->value(), 1.0);
}

TEST_F(AGroupTest, NestedListInitializationWithOneObject_Constructor){
    Group supergroup({Group({g1})});
    ASSERT_FALSE(supergroup.isLeaf());

    ASSERT_EQ(supergroup.representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0][0].representative()->value(), 1.1);
}

TEST_F(AGroupTest, NestedListInitialization_Constructor){
    Group supergroup({{g1, g2, g3}});
    ASSERT_FALSE(supergroup.isLeaf());

    ASSERT_EQ(supergroup.representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0][0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0][1].representative()->value(), 1.2);
    ASSERT_EQ(supergroup[0][2].representative()->value(), 1.0);
}

TEST_F(AGroupTest, ListInitialization_Empty){
    Group supergroup; // Careful, Group supergroup({}); will fail for intel compilers.
    ASSERT_EQ(supergroup.representative(), nullptr);
    ASSERT_TRUE(supergroup.isLeaf());
}

TEST_F(AGroupTest, EmplaceBack){
    Group supergroup;
    supergroup.emplace_back(g1);
    supergroup.emplace_back(g2);
    supergroup.emplace_back(g3);

    ASSERT_EQ(supergroup.representative()->value(), 1.1);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[1].representative()->value(), 1.2);
    ASSERT_EQ(supergroup[2].representative()->value(), 1.0);
}

TEST_F(AGroupTest, Sort){
    Group supergroup({g1, g2, g3});

    supergroup.sort();

    ASSERT_EQ(supergroup.representative()->value(), 1.0);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.0);
    ASSERT_EQ(supergroup[1].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[2].representative()->value(), 1.2);
}

TEST_F(AGroupTest, SortAll_Simple){
    Group supergroup({g1, g2, g3});

    supergroup.sortAll();

    ASSERT_EQ(supergroup.representative()->value(), 1.0);
    ASSERT_EQ(supergroup[0].representative()->value(), 1.0);
    ASSERT_EQ(supergroup[1].representative()->value(), 1.1);
    ASSERT_EQ(supergroup[2].representative()->value(), 1.2);
}

TEST_F(AGroupTest, SortAll_Nested){

    Group supergroup1({g2, g1});
    Group supergroup2({g2, g3});
    Group supersupergroup({supergroup1, supergroup2});

    supersupergroup.sortAll();

    ASSERT_EQ(supersupergroup.representative()->value(), 1.0);
    ASSERT_EQ(supersupergroup[0].representative()->value(), 1.0);
    ASSERT_EQ(supersupergroup[0][0].representative()->value(), 1.0);
    ASSERT_EQ(supersupergroup[0][1].representative()->value(), 1.2);

    ASSERT_EQ(supersupergroup[1].representative()->value(), 1.1);
    ASSERT_EQ(supersupergroup[1][0].representative()->value(), 1.1);
    ASSERT_EQ(supersupergroup[1][1].representative()->value(), 1.2);
}


TEST_F(AGroupTest, NumberOfLeaves) {
    Group supergroup1({g2, g1});
    Group supergroup2({g2, g3});
    Group supersupergroup({supergroup1, supergroup2});

    ASSERT_EQ(g1.numberOfLeaves(), 1);
    ASSERT_EQ(g2.numberOfLeaves(), 1);
    ASSERT_EQ(g3.numberOfLeaves(), 1);
    ASSERT_EQ(supergroup1.numberOfLeaves(), 2);
    ASSERT_EQ(supergroup2.numberOfLeaves(), 2);
    ASSERT_EQ(supersupergroup.numberOfLeaves(), 4);
}

TEST_F(AGroupTest, Merge) {
    Group supergroup0;
    ASSERT_EQ(supergroup0.representative(), nullptr);
    ASSERT_TRUE(supergroup0.isLeaf());

    Group supergroup1({g1});
    supergroup0 += supergroup1;//supergroup0.merge(supergroup1);

    ASSERT_EQ(supergroup0.representative()->value(), 1.1);
    ASSERT_EQ(supergroup0[0].representative()->value(), 1.1);

    Group supergroup2({g2,g3});
    supergroup0 += supergroup2;

    ASSERT_EQ(supergroup0.representative()->value(), 1.1);
    ASSERT_EQ(supergroup0[0].representative()->value(), 1.1);
    ASSERT_EQ(supergroup0[1].representative()->value(), 1.2);
    ASSERT_EQ(supergroup0[2].representative()->value(), 1.0);
}

TEST_F(AGroupTest, MakeSubgroup) {
    Group g({g1,g2,g3});

    ASSERT_EQ(g.size(),3);
    ASSERT_TRUE(g[0].isLeaf());
    ASSERT_TRUE(g[1].isLeaf());
    ASSERT_TRUE(g[2].isLeaf());
    ASSERT_EQ(g[0].representative()->value(),1.1);
    ASSERT_EQ(g[1].representative()->value(),1.2);
    ASSERT_EQ(g[2].representative()->value(),1.0);

    g.makeSubgroup({g.begin()+1, g.begin()}); // results in g: {g3,{g1,g2}}

    ASSERT_EQ(g.size(),2);
    ASSERT_TRUE(g[0].isLeaf());
    ASSERT_EQ(g[1].size(),2);
    ASSERT_EQ(g[0].representative()->value(),1.0);
    ASSERT_TRUE(g[1][0].isLeaf());
    ASSERT_EQ(g[1][0].representative()->value(),1.1);
    ASSERT_TRUE(g[1][1].isLeaf());
    ASSERT_EQ(g[1][1].representative()->value(),1.2);
}

TEST_F(AGroupTest, PlusOperator) {
    auto g = h1+h2;  // results in g: {g1, g2}}

    ASSERT_EQ(g.size(),2);
    ASSERT_TRUE(g[0].isLeaf());
    ASSERT_EQ(g[0].representative()->value(),1.1);
    ASSERT_TRUE(g[1].isLeaf());
    ASSERT_EQ(g[1].representative()->value(),1.2);

    g += h3; // results in g: {g1, g2, g3}
    ASSERT_EQ(g.size(),3);
    ASSERT_TRUE(g[0].isLeaf());
    ASSERT_EQ(g[0].representative()->value(),1.1);
    ASSERT_TRUE(g[1].isLeaf());
    ASSERT_EQ(g[1].representative()->value(),1.2);
    ASSERT_TRUE(g[2].isLeaf());
    ASSERT_EQ(g[2].representative()->value(),1.0);
}

TEST_F(AGroupTest, AllSampleIds) {
    Group g({g1,{g2,g3}});

    ASSERT_THAT(g.allSampleIds(), ElementsAre(0,1,2));

    ASSERT_THAT(g[0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(g[1].allSampleIds(), ElementsAre(1,2));

    ASSERT_THAT(g[1][0].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(g[1][1].allSampleIds(), ElementsAre(2));
}

TEST_F(AGroupTest, Print) {
    Group g({g1,{g2,g3}});
    std::stringstream ss;
    ss << g;

    ASSERT_EQ(ss.str(), std::string("{{0},{{1},{2}}}"));
}

TEST_F(AGroupTest, EmptyPrint) {
    Group g;
    std::stringstream ss;
    ss << g;
    ASSERT_EQ(ss.str(), std::string("{}"));
}

TEST_F(AGroupTest, Average) {
    Group g({g1,g3});

    auto averagedStructure = g.averagedMaximumPositionsVector();
    ASSERT_EQ(averagedStructure.weight, 2);
    Eigen::VectorXd avg = Eigen::VectorXd::Zero(
            TestMolecules::H2::ElectronsInCores::normal.electrons().positionsVector().numberOfEntities()
            * TestMolecules::H2::ElectronsInCores::normal.electrons().positionsVector().entityLength());

    ASSERT_TRUE(averagedStructure.positions.asEigenVector().isApprox(avg));
}

TEST_F(AGroupTest, Average_OneMember) {
    Group g({g1});

    auto averagedStructure = g.averagedMaximumPositionsVector();
    ASSERT_EQ(averagedStructure.weight, 1);

    ASSERT_TRUE(averagedStructure.positions.asEigenVector().isApprox(
            g1.representative()->maximum().positionsVector().asEigenVector()));
}