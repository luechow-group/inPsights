//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCGUI_STRINGMETHODVALUEPLOTTER_H
#define AMOLQCGUI_STRINGMETHODVALUEPLOTTER_H


#include <Qt3DCore/QEntity>
#include <QtCharts/QLineSeries>
#include "StringMethod.h"

class StringMethodValuesPlotter {
public:
    StringMethodValuesPlotter(){};

    QtCharts::QLineSeries* getLineSeries(const BSplines::ArcLengthParametrizedBSpline &arcLengthParametrizedBSpline,
                                         const unsigned resolution){

      QtCharts::QLineSeries *series = new QtCharts::QLineSeries();

      for (unsigned j = 0; j <= resolution; ++j) {
        double u = double(j) / double(resolution);
        auto p = arcLengthParametrizedBSpline.evaluate(u);
        series->append(u, p(0));
      }
      return series;
    }
};

#endif //AMOLQCGUI_STRINGMETHODVALUEPLOTTER_H
