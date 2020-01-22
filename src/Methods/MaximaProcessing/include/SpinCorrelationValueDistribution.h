/* Copyright 2020 Michael Heuer
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

#ifndef INPSIGHTS_SPINCORRELATIONVALUEDISTRIBUTION_H
#define INPSIGHTS_SPINCORRELATIONVALUEDISTRIBUTION_H

#include <Statistics.h>

class SpinCorrelationValueDistribution {
public:
    SpinCorrelationValueDistribution(Eigen::Index oneSidedNonzeroBinCount = 12);

    Eigen::Index calculateBinIndex(double spinCorrelation);

    void addSpinStatistic(const TriangularMatrixStatistics& spinCorrelations);

    Eigen::VectorXd getHistogramVector();

private:
    Eigen::Index oneSidedNonzeroBinCount_,binCount_;
    Eigen::VectorXd bins_;
};

#endif //INPSIGHTS_SPINCORRELATIONVALUEDISTRIBUTION_H
