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
#include <Cluster.h>
#include <Reference.h>
#include <TestMolecules.h>

using namespace testing;

class AClusterTest : public ::testing::Test {
public:
    Cluster g1,g2,g3, h1,h2,h3;
    AtomsVector atoms;

    void SetUp() override {
        atoms = TestMolecules::H2::ElectronsInCores::normal.atoms();
        g1 = Cluster({atoms, 1.1,TestMolecules::H2::ElectronsInCores::normal.electrons(), 0});
        g2 = Cluster({atoms, 1.2,TestMolecules::H2::ElectronsInCores::translated.electrons(),1});
        g3 = Cluster({atoms, 1.0,TestMolecules::H2::ElectronsInCores::flippedSpins.electrons(), 2});
        h1 = Cluster({g1});
        h2 = Cluster({g2});
        h3 = Cluster({g3});
    }
};


TEST_F(AClusterTest, DefaultConstructor) {
    Cluster supercluster0;
    ASSERT_EQ(supercluster0.representative(), nullptr);
    ASSERT_EQ(supercluster0.getSelectedElectronsCount(), 0);
    ASSERT_TRUE(supercluster0.isLeaf());
}

TEST_F(AClusterTest, CopyConstructor){
    Cluster cluster(g1);
    ASSERT_TRUE(cluster.isLeaf());
    ASSERT_EQ(cluster.representative()->value(), 1.1);
    ASSERT_EQ(cluster.getSelectedElectronsCount(), 2);
}

TEST_F(AClusterTest, Preallocation){
    Cluster supercluster(3);

    ASSERT_FALSE(supercluster.isLeaf());
    ASSERT_TRUE(supercluster[0].isLeaf());
    ASSERT_EQ(supercluster.getSelectedElectronsCount(), 0);

    ASSERT_EQ(supercluster.representative(), supercluster.front().representative());
    ASSERT_EQ(supercluster.front().representative(), nullptr);
    ASSERT_EQ(supercluster[1].representative(), nullptr);
    ASSERT_EQ(supercluster[2].representative(), nullptr);

    supercluster[0] = g1;
    supercluster[1] = g2;
    supercluster[2] = g3;

    ASSERT_FALSE(supercluster.isLeaf());
    ASSERT_EQ(supercluster.representative(), g1.representative());
    ASSERT_EQ(supercluster[0].representative(), g1.representative());
    ASSERT_EQ(supercluster.front().getSelectedElectronsCount(), 2);
}

TEST_F(AClusterTest, ListInitialization){
    Cluster supercluster = {g1, g2, g3};
    ASSERT_FALSE(supercluster.isLeaf());

    ASSERT_EQ(supercluster.representative()->value(), 1.1);
    ASSERT_EQ(supercluster.getSelectedElectronsCount(), 2);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[1].representative()->value(), 1.2);
    ASSERT_EQ(supercluster[2].representative()->value(), 1.0);
}

TEST_F(AClusterTest, ListInitialization_Constructor){
    Cluster supercluster({g1, g2, g3});
    ASSERT_FALSE(supercluster.isLeaf());

    ASSERT_EQ(supercluster.representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[1].representative()->value(), 1.2);
    ASSERT_EQ(supercluster[2].representative()->value(), 1.0);
}

TEST_F(AClusterTest, NestedListInitializationWithOneObject_Constructor){
    Cluster supercluster({Cluster({g1})});
    ASSERT_FALSE(supercluster.isLeaf());

    ASSERT_EQ(supercluster.representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0][0].representative()->value(), 1.1);
}

TEST_F(AClusterTest, NestedListInitialization_Constructor){
    Cluster supercluster({{g1, g2, g3}});
    ASSERT_FALSE(supercluster.isLeaf());

    ASSERT_EQ(supercluster.representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0][0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0][1].representative()->value(), 1.2);
    ASSERT_EQ(supercluster[0][2].representative()->value(), 1.0);
}

TEST_F(AClusterTest, ListInitialization_Empty){
    Cluster supercluster; // Careful, Cluster supercluster({}); will fail for intel compilers.
    ASSERT_EQ(supercluster.representative(), nullptr);
    ASSERT_TRUE(supercluster.isLeaf());
}

TEST_F(AClusterTest, EmplaceBack){
    Cluster supercluster;
    supercluster.emplace_back(g1);
    supercluster.emplace_back(g2);
    supercluster.emplace_back(g3);

    ASSERT_EQ(supercluster.representative()->value(), 1.1);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[1].representative()->value(), 1.2);
    ASSERT_EQ(supercluster[2].representative()->value(), 1.0);
}

TEST_F(AClusterTest, Sort){
    Cluster supercluster({g1, g2, g3});

    supercluster.sort();

    ASSERT_EQ(supercluster.representative()->value(), 1.0);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.0);
    ASSERT_EQ(supercluster[1].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[2].representative()->value(), 1.2);
}

TEST_F(AClusterTest, SortAll_Simple){
    Cluster supercluster({g1, g2, g3});

    supercluster.sortAll();

    ASSERT_EQ(supercluster.representative()->value(), 1.0);
    ASSERT_EQ(supercluster[0].representative()->value(), 1.0);
    ASSERT_EQ(supercluster[1].representative()->value(), 1.1);
    ASSERT_EQ(supercluster[2].representative()->value(), 1.2);
}

TEST_F(AClusterTest, SortAll_Nested){

    Cluster supercluster1({g2, g1});
    Cluster supercluster2({g2, g3});
    Cluster supersupercluster({supercluster1, supercluster2});

    supersupercluster.sortAll();

    ASSERT_EQ(supersupercluster.representative()->value(), 1.0);
    ASSERT_EQ(supersupercluster[0].representative()->value(), 1.0);
    ASSERT_EQ(supersupercluster[0][0].representative()->value(), 1.0);
    ASSERT_EQ(supersupercluster[0][1].representative()->value(), 1.2);

    ASSERT_EQ(supersupercluster[1].representative()->value(), 1.1);
    ASSERT_EQ(supersupercluster[1][0].representative()->value(), 1.1);
    ASSERT_EQ(supersupercluster[1][1].representative()->value(), 1.2);
}


TEST_F(AClusterTest, NumberOfLeaves) {
    Cluster supercluster1({g2, g1});
    Cluster supercluster2({g2, g3});
    Cluster supersupercluster({supercluster1, supercluster2});

    ASSERT_EQ(g1.numberOfLeaves(), 1);
    ASSERT_EQ(g2.numberOfLeaves(), 1);
    ASSERT_EQ(g3.numberOfLeaves(), 1);
    ASSERT_EQ(supercluster1.numberOfLeaves(), 2);
    ASSERT_EQ(supercluster2.numberOfLeaves(), 2);
    ASSERT_EQ(supersupercluster.numberOfLeaves(), 4);
}

TEST_F(AClusterTest, Merge) {
    Cluster supercluster0;
    ASSERT_EQ(supercluster0.representative(), nullptr);
    ASSERT_TRUE(supercluster0.isLeaf());

    Cluster supercluster1({g1});
    supercluster0 += supercluster1;//supercluster0.merge(supercluster1);

    ASSERT_EQ(supercluster0.representative()->value(), 1.1);
    ASSERT_EQ(supercluster0[0].representative()->value(), 1.1);

    Cluster supercluster2({g2,g3});
    supercluster0 += supercluster2;

    ASSERT_EQ(supercluster0.representative()->value(), 1.1);
    ASSERT_EQ(supercluster0[0].representative()->value(), 1.1);
    ASSERT_EQ(supercluster0[1].representative()->value(), 1.2);
    ASSERT_EQ(supercluster0[2].representative()->value(), 1.0);
}

TEST_F(AClusterTest, MakeSubcluster) {
    Cluster g({g1,g2,g3});

    ASSERT_EQ(g.size(),3);
    ASSERT_TRUE(g[0].isLeaf());
    ASSERT_TRUE(g[1].isLeaf());
    ASSERT_TRUE(g[2].isLeaf());
    ASSERT_EQ(g[0].representative()->value(),1.1);
    ASSERT_EQ(g[1].representative()->value(),1.2);
    ASSERT_EQ(g[2].representative()->value(),1.0);

    g.makeSubcluster({g.begin()+1, g.begin()}); // results in g: {g3,{g1,g2}}

    ASSERT_EQ(g.size(),2);
    ASSERT_TRUE(g[0].isLeaf());
    ASSERT_EQ(g[1].size(),2);
    ASSERT_EQ(g[0].representative()->value(),1.0);
    ASSERT_TRUE(g[1][0].isLeaf());
    ASSERT_EQ(g[1][0].representative()->value(),1.1);
    ASSERT_TRUE(g[1][1].isLeaf());
    ASSERT_EQ(g[1][1].representative()->value(),1.2);
}

TEST_F(AClusterTest, PlusOperator) {
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

TEST_F(AClusterTest, AllSampleIds) {
    Cluster g({g1,{g2,g3}});

    ASSERT_THAT(g.allSampleIds(), ElementsAre(0,1,2));

    ASSERT_THAT(g[0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(g[1].allSampleIds(), ElementsAre(1,2));

    ASSERT_THAT(g[1][0].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(g[1][1].allSampleIds(), ElementsAre(2));
}

TEST_F(AClusterTest, Print) {
    Cluster g({g1,{g2,g3}});
    std::stringstream ss;
    ss << g;

    ASSERT_EQ(ss.str(), std::string("{{0},{{1},{2}}}"));
}

TEST_F(AClusterTest, EmptyPrint) {
    Cluster g;
    std::stringstream ss;
    ss << g;
    ASSERT_EQ(ss.str(), std::string("{}"));
}