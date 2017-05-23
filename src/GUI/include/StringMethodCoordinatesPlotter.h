//
// Created by heuer on 12.05.17.
//

#ifndef AMOLQCGUI_BSPLINEPLOTTER_H
#define AMOLQCGUI_BSPLINEPLOTTER_H

#include "Polyline.h"
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

      for (int j = 0; j < pointsList.size(); ++j) {
        Polyline pl(root, Qt::red, pointsList[j], radius);
      }

    }


private:

};

#endif //AMOLQCGUI_BSPLINEPLOTTER_H
