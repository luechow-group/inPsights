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
#include <BestMatch.h>
#include <spdlog/spdlog.h>

using namespace testing;

TEST(ABestMatchTest, CombinePermutations){
    Eigen::PermutationMatrix<Eigen::Dynamic>p1,p2;
    p1.setIdentity(4);
    p2.setIdentity(4);

    Eigen::VectorXi expected(8);
    expected << 0,1,2,3,4,5,6,7;

    ASSERT_TRUE(p1.indices().base().isApprox(expected.segment(0,4)));
    ASSERT_TRUE(p2.indices().base().isApprox(expected.segment(0,4)));

    auto p = BestMatch::combinePermutations(p1,p2);
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

    auto p = BestMatch::combinePermutations(p1,p2);
    ASSERT_TRUE(p.indices().base().isApprox(expected));
}

TEST(ABestMatchTest, PermutationToFront) {
    std::list<long> indices;

    indices.emplace_back(3);
    indices.emplace_back(7);

    auto permutation = BestMatch::getPermutationToFront(indices,10);

    Eigen::VectorXi refIndices(10);
    refIndices << 2,3,4,0,5,6,7,1,8,9;

    ASSERT_EQ(permutation.indices(), refIndices);
};

TEST(ABestMatchTest, HeadToFullPermutation) {
    Eigen::VectorXi indices(3);
    indices << 2,0,1;

    Eigen::PermutationMatrix<Eigen::Dynamic> partialPermutation(indices);
    auto permutation = BestMatch::headToFullPermutation(partialPermutation, 10);

    Eigen::VectorXi refIndices(10);
    refIndices << 2,0,1,3,4,5,6,7,8,9;

    ASSERT_EQ(permutation.indices(), refIndices);
};

TEST(ABestMatchTest, TailToFullPermutation) {
    Eigen::VectorXi indices(3);
    indices << 2,0,1;

    Eigen::PermutationMatrix<Eigen::Dynamic> partialPermutation(indices);
    auto permutation = BestMatch::tailToFullPermutation(partialPermutation, 10);

    Eigen::VectorXi refIndices(10);
    refIndices << 0,1,2,3,4,5,6,9,7,8;

    ASSERT_EQ(permutation.indices(), refIndices);
};

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