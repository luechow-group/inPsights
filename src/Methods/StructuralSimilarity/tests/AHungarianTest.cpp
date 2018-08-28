//
// Created by Morian Sonnet on 22.05.2017.
//

#include <gmock/gmock.h>
#include "Hungarian.h"

TEST(HungarianTest, OneSolution) {
    Eigen::MatrixXd input(4,4);
    input << \
    88,83,69,92,\
    77,37,49,92,\
    11,69,5,86,\
    8,9,98,23;

    Eigen::MatrixXd expectedOutput(4,4);
    expectedOutput << \
    0,0,1,0,\
    0,1,0,0,\
    1,0,0,0,\
    0,0,0,1;

    Eigen::MatrixXd output(4,4);
    Hungarian::findMatching(input,output,MATCH_MIN);
    ASSERT_EQ(output,expectedOutput);
}


TEST(HungarianTest, TwoSolutions) {
    Eigen::MatrixXd input(4,4);
    input <<\
    90,75,75,80,\
    35,85,55,65,\
    125,95,90,105,\
    45,110,95,115;

    Eigen::MatrixXd output(4,4);
    Hungarian::findMatching(input,output,MATCH_MIN);

    Eigen::MatrixXd expectedOutputs[2] = {Eigen::MatrixXd(4,4),Eigen::MatrixXd(4,4)};
    expectedOutputs[0]<< \
    0,1,0,0,\
    0,0,0,1,\
    0,0,1,0,\
    1,0,0,0;

    expectedOutputs[1]<<\
    0,0,0,1,\
    0,0,1,0,\
    0,1,0,0,\
    1,0,0,0;

    ASSERT_THAT(output,::testing::AnyOf(expectedOutputs[0],expectedOutputs[1]));
}