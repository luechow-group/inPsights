// Copyright (C) 2018-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <PermutationHandling.h>
#include <spdlog/spdlog.h>

TEST(ABestMatchTest, CombinePermutations){
    Eigen::PermutationMatrix<Eigen::Dynamic>p1,p2;
    p1.setIdentity(4);
    p2.setIdentity(4);

    Eigen::VectorXi expected(8);
    expected << 0,1,2,3,4,5,6,7;

    ASSERT_TRUE(p1.indices().base().isApprox(expected.segment(0,4)));
    ASSERT_TRUE(p2.indices().base().isApprox(expected.segment(0,4)));

    auto p = PermutationHandling::combinePermutations(p1, p2);
    ASSERT_TRUE(p.indices().base().isApprox(expected));
}

TEST(ABestMatchTest, CombinePermutationsWithZeroLength){
    Eigen::PermutationMatrix<Eigen::Dynamic>p1,p2;
    p1.setIdentity(0);
    p2.setIdentity(4);

    Eigen::VectorXi expected(4);
    expected << 0,1,2,3;

    ASSERT_TRUE(p1.indices().base().isApprox(expected.segment(0,0)));
    ASSERT_TRUE(p2.indices().base().isApprox(expected.segment(0,4)));

    auto p = PermutationHandling::combinePermutations(p1, p2);
    ASSERT_TRUE(p.indices().base().isApprox(expected));
}

TEST(ABestMatchTest, PermutationToFront) {
    std::list<long> indices;

    indices.emplace_back(3);
    indices.emplace_back(7);

    auto permutation = PermutationHandling::getPermutationToFront(indices, 10);

    Eigen::VectorXi refIndices(10);
    refIndices << 2,3,4,0,5,6,7,1,8,9;

    ASSERT_EQ(permutation.indices(), refIndices);
};

TEST(ABestMatchTest, HeadToFullPermutation) {
    Eigen::VectorXi indices(3);
    indices << 2,0,1;

    Eigen::PermutationMatrix<Eigen::Dynamic> partialPermutation(indices);
    auto permutation = PermutationHandling::headToFullPermutation(partialPermutation, 10);

    Eigen::VectorXi refIndices(10);
    refIndices << 2,0,1,3,4,5,6,7,8,9;

    ASSERT_EQ(permutation.indices(), refIndices);
};

TEST(ABestMatchTest, TailToFullPermutation) {
    Eigen::VectorXi indices(3);
    indices << 2,0,1;

    Eigen::PermutationMatrix<Eigen::Dynamic> partialPermutation(indices);
    auto permutation = PermutationHandling::tailToFullPermutation(partialPermutation, 10);

    Eigen::VectorXi refIndices(10);
    refIndices << 0,1,2,3,4,5,6,9,7,8;

    ASSERT_EQ(permutation.indices(), refIndices);
};