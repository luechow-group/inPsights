//
// Created by Michael Heuer on 2019-02-02.
//
#include <gmock/gmock.h>
#include <GraphAnalysis.h>

using namespace testing;

class AGraphAnalysisTest : public ::testing::Test {
public:
    Eigen::MatrixXb A, B, C;

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
    };
};

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