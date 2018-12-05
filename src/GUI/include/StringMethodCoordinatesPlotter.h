//
// Created by heuer on 12.05.17.
//

#ifndef INPSIGHTS_BSPLINEPLOTTER_H
#define INPSIGHTS_BSPLINEPLOTTER_H

#include "Polyline.h"
#include "Sphere.h"
#include "ArcLengthParametrizedBSpline.h"

/* Plots a 3N dimensional spline by creating N 3D-dimensional spline
 * */

class StringMethodCoordinatesPlotter{
public:

    StringMethodCoordinatesPlotter(Qt3DCore::QEntity *root,
            const BSplines::ArcLengthParametrizedBSpline& arcLengthParametrizedBSpline,
                   const unsigned resolution, const float radius){

      long reducedDim = arcLengthParametrizedBSpline.getReducedDim();

      assert( reducedDim%3 == 0 && reducedDim > 0 );

      std::vector<std::vector<QVector3D>> pointsList(reducedDim/3);

      for (unsigned i = 0; i < resolution; ++i) {
            double u = double(i) / double(resolution - 1);
            Eigen::VectorXd evalResult = arcLengthParametrizedBSpline.reducedEvaluate(u, 0);

            for (int j = 0; j < pointsList.size(); ++j) {
                auto tmp = evalResult.segment(j*3,3);
                pointsList[j].push_back(QVector3D(float(tmp(0)),float(tmp(1)),float(tmp(2))));
            }
        }

      /* //markers for knots
      Eigen::VectorXd U = arcLengthParametrizedBSpline.getKnotVector();
      int p = arcLengthParametrizedBSpline.getDegree();
      int nCP = arcLengthParametrizedBSpline.getControlPointNumber();
      for (int j = 0+p; j < U.size()-p; ++j) {
        //Polyline pl(root, Spin::QColorFromType<Spin>(Spin::Alpha), pointsList[j], radius, true);
        Eigen::VectorXd uStructure= arcLengthParametrizedBSpline.reducedEvaluate(U(j), 0).tail(reducedDim);

        for (int k = 0; k < reducedDim/3; ++k) {
          QVector3D qVector3D;
          qVector3D.setX(uStructure(k*3+0));
          qVector3D.setY(uStructure(k*3+1));
          qVector3D.setZ(uStructure(k*3+2));
          new Sphere(root,Qt::green,qVector3D,0.03);
        }
      }*/

      for (int j = 0; j < pointsList.size()/2; ++j) {
        new Polyline(root, GuiHelper::QColorFromType<Spin>(Spin::alpha), pointsList[j], radius, true);
      }
      for (int j = pointsList.size()/2; j < pointsList.size(); ++j) {
       new Polyline(root, GuiHelper::QColorFromType<Spin>(Spin ::beta), pointsList[j], radius, true);
      }

    }


private:

};

#endif //INPSIGHTS_BSPLINEPLOTTER_H
