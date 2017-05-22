//
// Created by Moria on 22.05.2017.
//
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "Hungarian.h"
#include "Eigen/dense"
#include <iostream>

TEST(AHungarianTest, SimpleTests)
{
    Eigen::MatrixXd inputs[2]={Eigen::MatrixXd(4,4),Eigen::MatrixXd(4,4)};
    Eigen::MatrixXd outputsExpected[3]={Eigen::MatrixXd(4,4),Eigen::MatrixXd(4,4),Eigen::MatrixXd(4,4)};
    Eigen::MatrixXd outputs[2]={Eigen::MatrixXd(4,4),Eigen::MatrixXd(4,4)};
    inputs[0]<<88,83,69,92,77,37,49,92,11,69,5,86,8,9,98,23;
    inputs[1]<<90,75,75,80,35,85,55,65,125,95,90,105,45,110,95,115;
    outputsExpected[0]<<0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1;
    outputsExpected[1]<<0,1,0,0, 0,0,0,1, 0,0,1,0, 1,0,0,0;
    outputsExpected[2]<<0,0,0,1, 0,0,1,0, 0,1,0,0, 1,0,0,0;
    Hungarian::findMatching(inputs[0],outputs[0],MATCH_MIN);
    Hungarian::findMatching(inputs[1],outputs[1],MATCH_MIN);
    EXPECT_EQ(outputs[0],outputsExpected[0]);
    EXPECT_THAT(outputs[1],::testing::AnyOf(outputsExpected[1],outputsExpected[2]));
}
