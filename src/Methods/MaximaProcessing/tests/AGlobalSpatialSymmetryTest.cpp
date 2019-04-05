//
// Created by heuer on 05.04.19.
//

#include <gmock/gmock.h>
#include <GlobalSpatialSymmetrySorter.h>
#include <TestMolecules.h>

// needed?
#include <SimilarReferences.h>
#include <random>


using namespace testing;

class AGlobalClusterSorterTest : public ::testing::Test {
public:
    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        GlobalSimilaritySorter::settings.similarityRadius = 1; // prevent assert
    }
};