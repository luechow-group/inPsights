//
// Created by Michael Heuer on 30.06.16.
//

#include <gmock/gmock.h>
#include <Eigen/Core>

#include "BSplineStationaryPointFinder.h"
#include "BSpline.h"
#include "BSplineFromControlPolygon.h"
#include "ContainerConverter.h"

using namespace testing;

class ABSplineStationaryPointFinderTest : public Test {
public:
  void SetUp() override { }
};

TEST_F( ABSplineStationaryPointFinderTest, FindExtremaDegree3) {
  Eigen::Matrix<double,9,1> mat;
  mat << 0, 0, 0, 1, 0, 1, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs = myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsStationaryPointFinder(bs);

  std::vector<double> minima = bsStationaryPointFinder.getMinima(0);
  Eigen::VectorXd minimaEigen;
  minimaEigen.resize(minima.size());
  for (unsigned i = 0; i < minima.size(); ++i) { minimaEigen(i) = minima[i]; }

  std::vector<double> maxima = bsStationaryPointFinder.getMaxima(0);
  Eigen::VectorXd maximaEigen;
  maximaEigen.resize(maxima.size());
  for (unsigned i = 0; i < maxima.size(); ++i) { maximaEigen(i) = maxima[i]; }

  //ContainerConverter::printVectorXdForMathematica(minimaEigen); std::cout << std::endl;
  //ContainerConverter::printVectorXdForMathematica(maximaEigen); std::cout << std::endl;

  Eigen::VectorXd refMinima;
  refMinima.resize(3);
  refMinima << 0.0, 0.5, 1.0;
  ASSERT_TRUE(minimaEigen.isApprox(refMinima));

  Eigen::VectorXd refMaxima;
  refMaxima.resize(2);
  refMaxima << 1/3.0, 2/3.0;
  ASSERT_TRUE(maximaEigen.isApprox(refMaxima));
}

TEST_F( ABSplineStationaryPointFinderTest, FindMinimumInRepeatedRootDegree3) {
  Eigen::Matrix<double,9,1> mat;
  mat << 1, 1, 1, -1, -1, -1, 1, 1, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  std::vector<double> minima = bsExtremaFinder.getMinima(0);
  Eigen::VectorXd minimaEigen;
  minimaEigen.resize(minima.size());
  for (unsigned i = 0; i < minima.size(); ++i) { minimaEigen(i) = minima[i]; }

  std::vector<double> maxima = bsExtremaFinder.getMaxima(0);
  Eigen::VectorXd maximaEigen;
  maximaEigen.resize(maxima.size());
  for (unsigned i = 0; i < maxima.size(); ++i) { maximaEigen(i) = maxima[i]; }

  Eigen::VectorXd refMinima;
  refMinima.resize(1);
  refMinima << 0.5;
  ASSERT_TRUE(minimaEigen.isApprox(refMinima));

  Eigen::VectorXd refMaxima;
  refMaxima.resize(2);
  refMaxima << 0.0, 1.0;
  ASSERT_TRUE(maximaEigen.isApprox(refMaxima));
}

TEST_F( ABSplineStationaryPointFinderTest, FindMaximumInRepeatedRootDegree3) {
  Eigen::Matrix<double,9,1> mat;
  mat << -1,-1,-1, 1, 1, 1,-1,-1,-1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  std::vector<double> minima = bsExtremaFinder.getMinima(0);
  Eigen::VectorXd minimaEigen;
  minimaEigen.resize(minima.size());
  for (unsigned i = 0; i < minima.size(); ++i) { minimaEigen(i) = minima[i]; }

  std::vector<double> maxima = bsExtremaFinder.getMaxima(0);
  Eigen::VectorXd maximaEigen;
  maximaEigen.resize(maxima.size());
  for (unsigned i = 0; i < maxima.size(); ++i) { maximaEigen(i) = maxima[i]; }

  Eigen::VectorXd refMinima;
  refMinima.resize(2);
  refMinima << 0.0, 1.0;
  ASSERT_TRUE(minimaEigen.isApprox(refMinima));

  Eigen::VectorXd refMaxima;
  refMaxima.resize(1);
  refMaxima << 0.5;
  ASSERT_TRUE(maximaEigen.isApprox(refMaxima));
}

TEST_F( ABSplineStationaryPointFinderTest, FindBasinInRootIntervalDegree3) {
  Eigen::Matrix<double,6,1> mat;
  mat << 1,-1,-1,-1,-1, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  auto basins = bsExtremaFinder.getBasins(0);
  Eigen::MatrixXd basinsEigen;
  basinsEigen.resize(basins.size(),2);
  for (unsigned i = 0; i < basins.size(); ++i) {
    basinsEigen(i,0) = basins[i].first;
    basinsEigen(i,1) = basins[i].second;
  }

  auto plateaus = bsExtremaFinder.getPlateaus(0);
  Eigen::MatrixXd plateausEigen;
  plateausEigen.resize(plateaus.size(),2);
  for (unsigned i = 0; i < plateaus.size(); ++i) {
    plateausEigen(i,0) = plateaus[i].first;
    plateausEigen(i,1) = plateaus[i].second;
  }

  Eigen::MatrixXd refBasins;
  refBasins.resize(1,2);
  refBasins << 1/3.0, 2/3.0;
  ASSERT_TRUE(basinsEigen.isApprox(refBasins));

  Eigen::MatrixXd refPlateaus;
  refPlateaus.resize(0,2);
  ASSERT_TRUE(plateausEigen.isApprox(refPlateaus));
}

TEST_F( ABSplineStationaryPointFinderTest, FindPlateauInRootIntervalDegree3) {
  Eigen::Matrix<double,6,1> mat;
  mat << -1, 1, 1, 1, 1,-1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  auto basins = bsExtremaFinder.getBasins(0);
  Eigen::MatrixXd basinsEigen;
  basinsEigen.resize(basins.size(),2);
  for (unsigned i = 0; i < basins.size(); ++i) {
    basinsEigen(i,0) = basins[i].first;
    basinsEigen(i,1) = basins[i].second;
  }

  auto plateaus = bsExtremaFinder.getPlateaus(0);
  Eigen::MatrixXd plateausEigen;
  plateausEigen.resize(plateaus.size(),2);
  for (unsigned i = 0; i < plateaus.size(); ++i) {
    plateausEigen(i,0) = plateaus[i].first;
    plateausEigen(i,1) = plateaus[i].second;
  }

  Eigen::MatrixXd refBasins;
  refBasins.resize(0,2);
  ASSERT_TRUE(basinsEigen.isApprox(refBasins));

  Eigen::MatrixXd refPlateaus;
  refPlateaus.resize(1,2);
  refPlateaus << 1/3.0, 2/3.0;
  ASSERT_TRUE(plateausEigen.isApprox(refPlateaus));
}

TEST_F( ABSplineStationaryPointFinderTest, FindReorientationInRootIntervalDegree3) {
  Eigen::Matrix<double,6,1> mat;
  mat << -1,-1,-1,-1,-1,-1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  auto basins = bsExtremaFinder.getBasins(0);
  Eigen::MatrixXd basinsEigen;
  basinsEigen.resize(basins.size(),2);
  for (unsigned i = 0; i < basins.size(); ++i) {
    basinsEigen(i,0) = basins[i].first;
    basinsEigen(i,1) = basins[i].second;
  }

  auto plateaus = bsExtremaFinder.getPlateaus(0);
  Eigen::MatrixXd plateausEigen;
  plateausEigen.resize(plateaus.size(),2);
  for (unsigned i = 0; i < plateaus.size(); ++i) {
    plateausEigen(i,0) = plateaus[i].first;
    plateausEigen(i,1) = plateaus[i].second;
  }

  auto reorientations = bsExtremaFinder.getPureReorientation(0);
  Eigen::MatrixXd reorientationsEigen;
  reorientationsEigen.resize(reorientations.size(),2);
  for (unsigned i = 0; i < reorientations.size(); ++i) {
    reorientationsEigen(i,0) = reorientations[i].first;
    reorientationsEigen(i,1) = reorientations[i].second;
  }

  Eigen::MatrixXd refBasins;
  refBasins.resize(0,2);
  ASSERT_TRUE(basinsEigen.isApprox(refBasins));

  Eigen::MatrixXd refPlateaus;
  refPlateaus.resize(0,2);
  ASSERT_TRUE(plateausEigen.isApprox(refPlateaus));

  Eigen::MatrixXd refReorientations;
  refReorientations.resize(1,2);
  refReorientations << 0.0, 1.0;
  ASSERT_TRUE(reorientationsEigen.isApprox(refReorientations));
}

TEST_F( ABSplineStationaryPointFinderTest, FindBasinInRootIntervalAtHeadAndTailDegree3) {
  Eigen::Matrix<double,9,1> mat;
  mat << -1,-1,-1,-1, 1,-1,-1,-1,-1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  std::vector<double> minima = bsExtremaFinder.getMinima(0);
  Eigen::VectorXd minimaEigen;
  minimaEigen.resize(minima.size());
  for (unsigned i = 0; i < minima.size(); ++i) { minimaEigen(i) = minima[i]; }

  std::vector<double> maxima = bsExtremaFinder.getMaxima(0);
  Eigen::VectorXd maximaEigen;
  maximaEigen.resize(maxima.size());
  for (unsigned i = 0; i < maxima.size(); ++i) { maximaEigen(i) = maxima[i]; }


  auto basins = bsExtremaFinder.getBasins(0);
  Eigen::MatrixXd basinsEigen;
  basinsEigen.resize(basins.size(),2);
  for (unsigned i = 0; i < basins.size(); ++i) {
    basinsEigen(i,0) = basins[i].first;
    basinsEigen(i,1) = basins[i].second;
  }

  auto plateaus = bsExtremaFinder.getPlateaus(0);
  Eigen::MatrixXd plateausEigen;
  plateausEigen.resize(plateaus.size(),2);
  for (unsigned i = 0; i < plateaus.size(); ++i) {
    plateausEigen(i,0) = plateaus[i].first;
    plateausEigen(i,1) = plateaus[i].second;
  }

  auto reorientations = bsExtremaFinder.getPureReorientation(0);
  Eigen::MatrixXd reorientationsEigen;
  reorientationsEigen.resize(reorientations.size(),2);
  for (unsigned i = 0; i < reorientations.size(); ++i) {
    reorientationsEigen(i,0) = reorientations[i].first;
    reorientationsEigen(i,1) = reorientations[i].second;
  }

  Eigen::VectorXd refMinima;
  refMinima.resize(0);
  ASSERT_TRUE(minimaEigen.isApprox(refMinima));

  Eigen::VectorXd refMaxima;
  refMaxima.resize(1);
  refMaxima << 0.5;
  ASSERT_TRUE(maximaEigen.isApprox(refMaxima));

  Eigen::MatrixXd refBasins;
  refBasins.resize(2,2);
  refBasins << 0.0, 1/6.0, 5/6.0, 1.0;
  ASSERT_TRUE(basinsEigen.isApprox(refBasins));

  Eigen::MatrixXd refPlateaus;
  refPlateaus.resize(0,2);
  ASSERT_TRUE(plateausEigen.isApprox(refPlateaus));

  Eigen::MatrixXd refReorientations;
  refReorientations.resize(2,2);
  refReorientations << 0.0, 1/6.0, 5/6.0, 1.0;
  ASSERT_TRUE(reorientationsEigen.isApprox(refReorientations));
}

TEST_F( ABSplineStationaryPointFinderTest, FindSaddlePointDegree3) {
  Eigen::Matrix<double,5,1> mat;
  mat << 0, 1, 1, 1, 2;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();
  BSplineStationaryPointFinder bsExtremaFinder(bs);

  std::vector<double> minima = bsExtremaFinder.getMinima(0);
  Eigen::VectorXd minimaEigen;
  minimaEigen.resize(minima.size());
  for (unsigned i = 0; i < minima.size(); ++i) {
    std::cout << minima[i] << std::endl;
    minimaEigen(i) = minima[i]; }

  std::vector<double> maxima = bsExtremaFinder.getMaxima(0);
  Eigen::VectorXd maximaEigen;
  maximaEigen.resize(maxima.size());
  for (unsigned i = 0; i < maxima.size(); ++i) { maximaEigen(i) = maxima[i]; }

  std::vector<double> saddlePoints = bsExtremaFinder.getIncSaddlePoints(0);
  Eigen::VectorXd saddlePointsEigen;
  saddlePointsEigen.resize(saddlePoints.size());
  for (unsigned i = 0; i < saddlePoints.size(); ++i) { saddlePointsEigen(i) = saddlePoints[i]; }


  auto basins = bsExtremaFinder.getBasins(0);
  Eigen::MatrixXd basinsEigen;
  basinsEigen.resize(basins.size(),2);
  for (unsigned i = 0; i < basins.size(); ++i) {
    basinsEigen(i,0) = basins[i].first;
    basinsEigen(i,1) = basins[i].second;
  }

  auto plateaus = bsExtremaFinder.getPlateaus(0);
  Eigen::MatrixXd plateausEigen;
  plateausEigen.resize(plateaus.size(),2);
  for (unsigned i = 0; i < plateaus.size(); ++i) {
    plateausEigen(i,0) = plateaus[i].first;
    plateausEigen(i,1) = plateaus[i].second;
  }

  auto reorientations = bsExtremaFinder.getPureReorientation(0);
  Eigen::MatrixXd reorientationsEigen;
  reorientationsEigen.resize(reorientations.size(),2);
  for (unsigned i = 0; i < reorientations.size(); ++i) {
    reorientationsEigen(i,0) = reorientations[i].first;
    reorientationsEigen(i,1) = reorientations[i].second;
  }

  Eigen::VectorXd refMinima;
  refMinima.resize(0);
  ASSERT_TRUE(minimaEigen.isApprox(refMinima));

  Eigen::VectorXd refMaxima;
  refMaxima.resize(0);
  ASSERT_TRUE(maximaEigen.isApprox(refMaxima));

  Eigen::VectorXd refSaddlePoints;
  refSaddlePoints.resize(1);
  refSaddlePoints << 0.5;
  ASSERT_TRUE(saddlePointsEigen.isApprox(refSaddlePoints));

  Eigen::MatrixXd refBasins;
  refBasins.resize(0,2);
  ASSERT_TRUE(basinsEigen.isApprox(refBasins));

  Eigen::MatrixXd refPlateaus;
  refPlateaus.resize(0,2);
  ASSERT_TRUE(plateausEigen.isApprox(refPlateaus));

  Eigen::MatrixXd refReorientations;
  refReorientations.resize(0,2);
  ASSERT_TRUE(reorientationsEigen.isApprox(refReorientations));
}