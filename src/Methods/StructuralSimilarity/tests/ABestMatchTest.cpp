// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <BestMatch.h>

TEST(ABestMatchTest, LesserOperator){
    Eigen::VectorXi indices(3);
    indices << 0,1,2;
    Eigen::PermutationMatrix<Eigen::Dynamic> perm(indices);

    BestMatch::AscendingMetricResult a = {0,perm};
    BestMatch::AscendingMetricResult b = {1,perm};
    ASSERT_TRUE(a<b);

    BestMatch::DescendingMetricResult c= {0,perm};
    BestMatch::DescendingMetricResult  d= {1,perm};

    ASSERT_TRUE(d<c);
}

TEST(ABestMatchTest, Sort){
    Eigen::VectorXi indices(3);
    indices << 0,1,2;
    Eigen::PermutationMatrix<Eigen::Dynamic> perm(indices);

    BestMatch::AscendingMetricResult a = {0,perm};
    BestMatch::AscendingMetricResult b = {1,perm};
    BestMatch::AscendingMetricResult c = {2,perm};

    std::vector< BestMatch::AscendingMetricResult> vec = {a,c,b};
    std::sort(vec.begin(),vec.end());

}

TEST(ABestMatchTest, ComplexAscendingMetricSort) {

    Eigen::VectorXi p1(3), p2(3), p3(3), p4(3), p5(3), p6(3);
    p1 << 0,1,2;
    p2 << 0,2,1;
    p3 << 1,2,0;
    p4 << 1,0,2;
    p5 << 2,0,1;
    p6 << 2,1,0;

    std::vector<BestMatch::AscendingMetricResult > results = {
            {1.0, Eigen::PermutationMatrix<Eigen::Dynamic>(p2)},
            {1.0, Eigen::PermutationMatrix<Eigen::Dynamic>(p1)},
            {0.9, Eigen::PermutationMatrix<Eigen::Dynamic>(p5)},
            {0.8, Eigen::PermutationMatrix<Eigen::Dynamic>(p4)},
            {0.9, Eigen::PermutationMatrix<Eigen::Dynamic>(p3)},
            {1.0, Eigen::PermutationMatrix<Eigen::Dynamic>(p6)},
    };
    std::sort(results.begin(), results.end());

    std::vector<Eigen::VectorXi> expectedPermIndices(6, Eigen::VectorXi(3));
    expectedPermIndices[0] = p4;
    expectedPermIndices[1] = p3;
    expectedPermIndices[2] = p5;
    expectedPermIndices[3] = p1;
    expectedPermIndices[4] = p2;
    expectedPermIndices[5] = p6;

    ASSERT_EQ(results.size(), expectedPermIndices.size());
    for (size_t i = 0; i < results.size(); ++i) {
        ASSERT_EQ(results[i].permutation.indices(), expectedPermIndices[i]);
    }
}

TEST(ABestMatchTest, ComplexDescendingMetricSort) {
    Eigen::VectorXi p1(3), p2(3), p3(3), p4(3), p5(3), p6(3);
    p1 << 0,1,2;
    p2 << 0,2,1;
    p3 << 1,2,0;
    p4 << 1,0,2;
    p5 << 2,0,1;
    p6 << 2,1,0;


    std::vector<BestMatch::DescendingMetricResult > results = {
            {1.0, Eigen::PermutationMatrix<Eigen::Dynamic>(p2)},
            {1.0, Eigen::PermutationMatrix<Eigen::Dynamic>(p1)},
            {0.9, Eigen::PermutationMatrix<Eigen::Dynamic>(p5)},
            {0.8, Eigen::PermutationMatrix<Eigen::Dynamic>(p4)},
            {0.9, Eigen::PermutationMatrix<Eigen::Dynamic>(p3)},
            {1.0, Eigen::PermutationMatrix<Eigen::Dynamic>(p6)},
    };
    std::sort(results.begin(), results.end());

    std::vector<Eigen::VectorXi> expectedPermIndices(6, Eigen::VectorXi(3));
    expectedPermIndices[0] = p1;
    expectedPermIndices[1] = p2;
    expectedPermIndices[2] = p6;
    expectedPermIndices[3] = p3;
    expectedPermIndices[4] = p5;
    expectedPermIndices[5] = p4;

    ASSERT_EQ(results.size(), expectedPermIndices.size());
    for (size_t i = 0; i < results.size(); ++i) {
        ASSERT_EQ(results[i].permutation.indices(), expectedPermIndices[i]);
    }
}