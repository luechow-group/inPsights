//
// Created by Michael Heuer on 04.09.18.
//

#include "Comparators.h"
#include <TestMolecules.h>
#include <gtest/gtest.h>

TEST(AComparatorTest, Comparator) {

    auto ref1 = Reference(TestMolecules::threeElectrons::normal.electrons(),1.0);
    auto ref2 = ref1;
    ref2.maximum_.positionsVector().entity(0).translate({0,0, 0.1});
    ref2.maximum_.positionsVector().entity(1).translate({0,0,-0.1});

    Eigen::VectorXi perm(2);
    perm << 1,0;

    ref2.maximum_.slice({0,2}).permute(Eigen::PermutationMatrix<Eigen::Dynamic>(perm));
    ref2.negLogSqrdProbabilityDensity_ += 0.1;
    // deltaDistance = 0.1, deltaValue = 0.1;

    auto eps = std::numeric_limits<double>::epsilon();
    ASSERT_TRUE(Comparators::ValueEuclideanDistanceComparator(0.30,0.20)(ref1,ref2));
    ASSERT_TRUE(Comparators::ValueEuclideanDistanceComparator(0.10+eps,0.10+eps)(ref1,ref2));

    ASSERT_FALSE(Comparators::ValueEuclideanDistanceComparator(0.09,0.20)(ref1,ref2));
    ASSERT_FALSE(Comparators::ValueEuclideanDistanceComparator(0.20,0.09)(ref1,ref2));
    ASSERT_FALSE(Comparators::ValueEuclideanDistanceComparator(0.09,0.09)(ref1,ref2));
}
