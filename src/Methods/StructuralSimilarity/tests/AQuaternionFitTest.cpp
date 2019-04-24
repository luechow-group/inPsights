/* Copyright (C) 2017-2019 Michael Heuer.
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
#include <Eigen/Core>
#include <QuaternionFit.h>

using namespace testing;
using namespace Eigen;

TEST(AQuaternionFitTest, IdenticalGeometries) {
  MatrixXd refMat(4,3);
  MatrixXd fitMat(4,3);

  refMat << \
  1,0,0,\
  0,1,0,\
  0,0,1,\
  1,1,1;

  fitMat << \
  1,0,0,\
  0,1,0,\
  0,0,1,\
  1,1,1;

  QuaternionFit quatFit = QuaternionFit(refMat,fitMat);
  ASSERT_TRUE( quatFit.getRotationMatrix().isApprox(Eigen::Matrix3d::Identity()) );
}

TEST(AQuaternionFitTest, DISABLED_ThreeCollinearPoints) {
  ASSERT_TRUE(false); //TODO
}

TEST(AQuaternionFitTest, DISABLED_CoplanarPoints) {
  ASSERT_TRUE(false); //TODO
}

TEST(AQuaternionFitTest, DISABLED_TwoPoints) {
  ASSERT_TRUE(false); //TODO
}

TEST(AQuaternionFitTest, DISABLED_SinglePoint) {
  ASSERT_TRUE(false); //TODO
}

TEST(AQuaternionFitTest, RotateIdenticalGeometriesBy90DegreeAroundZAxis) {
  MatrixXd refMat(4, 3);
  MatrixXd fitMat(4, 3);

  refMat <<
    1, 0, 0,\
    0, 1, 0,\
    0, 0, 1,\
    1, 0, 1;

  fitMat <<
     0, 1, 0,\
    -1, 0, 0,\
     0, 0, 1,\
     0, 1, 1;

  QuaternionFit quatFit = QuaternionFit(refMat, fitMat);
  Matrix3d refRotMat;
  refRotMat << \
   0,1,0,\
  -1,0,0,\
   0,0,1;

  ASSERT_TRUE(quatFit.getRotationMatrix().isApprox(refRotMat));
}

TEST(AQuaternionFitTest, DifferentGeometriesRMSD) {
  // Acetylacetone geometry

  MatrixXd refMat(15,3);
  MatrixXd fitMat(15,3);

  refMat <<  0.837039463013,  0.016788427320, -2.525595956526,
             1.543513550153, -0.813089574597, -2.549514746791,
             1.414756694727,  0.939615423053, -2.588330890913,
             0.157123350877, -0.049408205921, -3.369553256885,
             0.060071572521, -0.006488221348, -1.238960132765,
             0.794033293747,  0.018022065819,  0.000474852563,
             1.872442556461,  0.045447119532, -0.009796692884,
             0.134542340451,  0.008553783896,  1.197861199576,
             0.807885933698,  0.034713426007,  2.524751492477,
             1.888410774041,  0.057732055006,  2.416018936208,
             0.478749565998,  0.911663008528,  3.081848148239,
             0.517505535617, -0.846061511818,  3.097060472785,
            -1.188496502471, -0.042294663829, -1.261659669635,
            -1.190036824170, -0.024442923129,  1.268285673962,
            -1.499835119606, -0.039635177950,  0.311246144880;

  fitMat << -2.5344764741,   -0.7559439560,    0.0000783747,
            -2.5850554288,   -1.3994142029,    0.8788953171,
            -2.5852111247,   -1.3990212850,   -0.8790195729,
            -3.3707668752,   -0.0634195032,    0.0002860078,
            -1.2393081666,    0.0071159529,    0.0000891199,
            -0.0069259062,   -0.7393224246,    0.0001171426,
            -0.0270888432,   -1.8179608046,    0.0001766557,
             1.1963512412,   -0.0903665930,    0.0000661919,
             2.5179086002,   -0.7746250207,    0.0000823538,
             2.4011372771,   -1.8545542125,    0.0001544348,
             3.0848864988,   -0.4690485428,   -0.8790550581,
             3.0849329122,   -0.4689320016,    0.8791495949,
            -1.2495714780,    1.2563573297,    0.0000552841,
             1.2782591621,    1.2338842886,   -0.0000066779,
             0.3235786053,    1.5521309757,    0.0000008317;

  QuaternionFit quatFit = QuaternionFit(refMat,fitMat);
}