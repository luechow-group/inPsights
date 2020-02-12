/* Copyright (C) 2019 Michael Heuer.
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
#include <GraphAnalysis.h>

using namespace testing;

class AGraphAnalysisTest : public ::testing::Test {
public:
    Eigen::MatrixXb A, B, C;
    Eigen::MatrixXd distMat;

    void SetUp() override {
        A = Eigen::MatrixXb(4, 4);
        // Two pairs: (0-1,) (2-3)
        A << 0, 1, 0, 0, \
             1, 0, 0, 0, \
             0, 0, 0, 1, \
             0, 0, 1, 0;

        // Cyclic and isolated: (0-1, 0-3, 1-3) (2)
        B = Eigen::MatrixXb(4, 4);
        B << 0, 1, 0, 1, \
             1, 0, 0, 1, \
             0, 0, 0, 0, \
             1, 1, 0, 0;

        // Isolated and selfconnected: (0) (1) (2) (3)
        C = Eigen::MatrixXb(4, 4);
        C << 0, 0, 0, 0, \
             0, 0, 0, 0, \
             0, 0, 1, 0, \
             0, 0, 0, 1;

        distMat = Eigen::MatrixXd(4, 4);
        distMat << 0.0, 0.0, 0.2, 0.4, \
                   0.0, 0.0, 0.4, 0.8, \
                   0.2, 0.4, 0.0, 1.2, \
                   0.4, 0.8, 1.2, 0.0;
    };
};

TEST_F(AGraphAnalysisTest, FindConnectedVertices) {
    Eigen::MatrixXb expected02 (4,4);
    expected02 << 1, 1, 1, 0, \
                  1, 1, 0, 0, \
                  1, 0, 1, 0, \
                  0, 0, 0, 1;
    ASSERT_TRUE(GraphAnalysis::lowerOrEqualFilter(distMat, 0.2).isApprox(expected02));

    Eigen::MatrixXb expected04 (4,4);
    expected04 << 1, 1, 1, 1, \
                  1, 1, 1, 0, \
                  1, 1, 1, 0, \
                  1, 0, 0, 1;
    ASSERT_TRUE(GraphAnalysis::lowerOrEqualFilter(distMat, 0.4).isApprox(expected04));

    Eigen::MatrixXb expected08 (4,4);
    expected08 << 1, 1, 1, 1, \
                  1, 1, 1, 1, \
                  1, 1, 1, 0, \
                  1, 1, 0, 1;
    ASSERT_TRUE( GraphAnalysis::lowerOrEqualFilter(distMat, 0.8).isApprox(expected08));

    Eigen::MatrixXb expected12 (4,4);
    expected12 << 1, 1, 1, 1, \
                  1, 1, 1, 1, \
                  1, 1, 1, 1, \
                  1, 1, 1, 1;
    ASSERT_TRUE(GraphAnalysis::lowerOrEqualFilter(distMat, 1.2).isApprox(expected12));
}

TEST_F(AGraphAnalysisTest, LowerOrEqualFilter) {
    auto adjacencyMat02 = GraphAnalysis::lowerOrEqualFilter(distMat, 0.2);
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat02, 0), ElementsAre(0, 1, 2));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat02, 1), ElementsAre(0, 1, 2));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat02, 2), ElementsAre(0, 1, 2));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat02, 3), ElementsAre(3));

    auto adjacencyMat04 = GraphAnalysis::lowerOrEqualFilter(distMat, 0.4);
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat04, 0), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat04, 1), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat04, 2), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat04, 3), ElementsAre(0, 1, 2, 3));

    auto adjacencyMat08 = GraphAnalysis::lowerOrEqualFilter(distMat, 0.8);
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat08, 0), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat08, 1), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat08, 2), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat08, 3), ElementsAre(0, 1, 2, 3));

    auto adjacencyMat12 = GraphAnalysis::lowerOrEqualFilter(distMat, 1.2);
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat12, 0), ElementsAre(0,1,2,3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat12, 1), ElementsAre(0,1,2,3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat12, 2), ElementsAre(0,1,2,3));
    ASSERT_THAT(GraphAnalysis::findConnectedVertices(adjacencyMat12, 3), ElementsAre(0,1,2,3));
}

TEST_F(AGraphAnalysisTest, TwoPairs) {
    auto clusters = GraphAnalysis::findGraphClusters(A);
    ASSERT_THAT(clusters, SizeIs(2));
    ASSERT_THAT(clusters[0], ElementsAre(0, 1));
    ASSERT_THAT(clusters[1], ElementsAre(2, 3));
}

TEST_F(AGraphAnalysisTest, CyclicAndIsolated) {
    auto clusters = GraphAnalysis::findGraphClusters(B);
    ASSERT_THAT(clusters, SizeIs(2));
    ASSERT_THAT(clusters[0], ElementsAre(0, 1, 3));
    ASSERT_THAT(clusters[1], ElementsAre(2));
}

TEST_F(AGraphAnalysisTest, IsolatedAndSelfconnected) {
    auto clusters = GraphAnalysis::findGraphClusters(C);
    ASSERT_THAT(clusters, SizeIs(4));
    ASSERT_THAT(clusters[0], ElementsAre(0));
    ASSERT_THAT(clusters[1], ElementsAre(1));
    ASSERT_THAT(clusters[2], ElementsAre(2));
    ASSERT_THAT(clusters[3], ElementsAre(3));
}

TEST_F(AGraphAnalysisTest, Filter) {
    Eigen::MatrixXd m(2,2);
    m << 0,0.45,0.5,1;

    Eigen::MatrixXb expected(2,2);
    expected << 0,0,1,1;

    ASSERT_EQ(GraphAnalysis::filter(m,0.5), expected);

    EXPECT_DEATH(GraphAnalysis::filter(m,1.1),"");
    m = m.array() - 1.0;
    EXPECT_DEATH(GraphAnalysis::filter(m,0.5),"");
}


TEST_F(AGraphAnalysisTest, FindSubsets) {

    std::vector<std::list<Eigen::Index>> subsets{{0, 1}, {2}, {3}};
    std::vector<std::list<Eigen::Index>> referenceSets{{0, 1, 2}, {3}};

    auto map = GraphAnalysis::findMergeMap(subsets, referenceSets);

    ASSERT_THAT(map,Contains(std::pair<Eigen::Index, Eigen::Index>(0,0)));
    ASSERT_THAT(map,Contains(std::pair<Eigen::Index, Eigen::Index>(1,0)));
    ASSERT_THAT(map,Contains(std::pair<Eigen::Index, Eigen::Index>(2,1)));
}

TEST_F(AGraphAnalysisTest, FindSubsetsUnorderd) {

    std::vector<std::list<Eigen::Index>> subsets{{2}, {3}, {0, 1}};
    std::vector<std::list<Eigen::Index>> referenceSets{{3}, {0, 1, 2}};

    auto map = GraphAnalysis::findMergeMap(subsets, referenceSets);

    ASSERT_THAT(map,Contains(std::pair<Eigen::Index, Eigen::Index>(0,1)));
    ASSERT_THAT(map,Contains(std::pair<Eigen::Index, Eigen::Index>(1,0)));
    ASSERT_THAT(map,Contains(std::pair<Eigen::Index, Eigen::Index>(2,1)));
}

TEST_F(AGraphAnalysisTest, FindKeysToValueInMap) {

    std::map<int, std::string> map = {
            {0, "b"},
            {1, "a"},
            {2, "b"}};

    std::vector<int> keysMappedToValueA, keysMappedToValueB, keysMappedToValueC;

    auto foundA = GraphAnalysis::findByValue(keysMappedToValueA, map, std::string("a"));
    auto foundB = GraphAnalysis::findByValue(keysMappedToValueB, map, std::string("b"));
    auto foundC = GraphAnalysis::findByValue(keysMappedToValueC, map, std::string("c"));

    ASSERT_TRUE(foundA);
    ASSERT_TRUE(foundB);
    ASSERT_FALSE(foundC);

    ASSERT_THAT(keysMappedToValueA,ElementsAre(1));
    ASSERT_THAT(keysMappedToValueB,ElementsAre(0,2));
    ASSERT_THAT(keysMappedToValueC,IsEmpty());
}

TEST_F(AGraphAnalysisTest, ReferenceSetWithoutMatchingSubsetDeath) {
    std::vector<std::list<Eigen::Index>> subsets{{0,1}, {2}, {3}};
    std::vector<std::list<Eigen::Index>> referenceSets{{0,1,2},{3},{4}};

    EXPECT_DEATH(GraphAnalysis::findMergeMap(subsets, referenceSets), "");
}

TEST_F(AGraphAnalysisTest, SubsetNotFoundInReferenceSetDeath) {
    std::vector<std::list<Eigen::Index>> subsets{{0,1}, {2}, {3}, {4}};
    std::vector<std::list<Eigen::Index>> referenceSets{{0,1,2},{3}};

    EXPECT_DEATH(GraphAnalysis::findMergeMap(subsets, referenceSets), "");
}
