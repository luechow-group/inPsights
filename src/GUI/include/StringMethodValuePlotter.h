//
// Created by heuer on 23.05.17.
//

#ifndef INPSIGHTS_STRINGMETHODVALUEPLOTTER_H
#define INPSIGHTS_STRINGMETHODVALUEPLOTTER_H


#include <Qt3DCore/QEntity>
#include <QtCharts/QLineSeries>
#include "StringMethod.h"

class StringMethodValuesPlotter {
public:
    StringMethodValuesPlotter() {};

    QtCharts::QLineSeries *getLineSeries(const BSplines::ArcLengthParametrizedBSpline &arcLengthParametrizedBSpline,
                                         const unsigned resolution) {

        QtCharts::QLineSeries *series = new QtCharts::QLineSeries();

        for (unsigned j = 0; j <= resolution; ++j) {
            double u = double(j) / double(resolution);
            auto p = arcLengthParametrizedBSpline.evaluate(u);
            series->append(u, p(0));
        }
        return series;
    }
};

#endif //INPSIGHTS_STRINGMETHODVALUEPLOTTER_H
