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
