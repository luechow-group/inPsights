/* Copyright (C) 2018-2019 Michael Heuer.
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

#include "Sinkhorn.h"

using namespace SOAP;

Eigen::MatrixXd Sinkhorn::Pgamma(const Eigen::MatrixXd &C, double gamma, double eps) {
    auto N = C.rows();
    auto M = C.cols();

    Eigen::VectorXd u = Eigen::VectorXd::Constant(N, 1. / double(N));
    Eigen::VectorXd v = Eigen::VectorXd::Constant(M, 1. / double(M));
    Eigen::VectorXd uOld(N), vOld(M);
    Eigen::MatrixXd K(N, M);

    double aN = 1.0 / double(N), aM = 1.0 / double(M);
    double uDeviation, vDeviation, err = std::pow(eps, 2);

    double lambda = 1. / gamma;
    K = (lambda * (C - Eigen::MatrixXd::Ones(N, M))).array().exp();

    do {
        uOld = u;
        uDeviation = 0.0;

        for (Eigen::Index i = 0; i < N; ++i) {
            u[i] = K.row(i) * v;
            uDeviation += pow(aN - uOld[i] * u[i], 2);
        }
        u = aN * u.array().cwiseInverse();

        vOld = v;
        vDeviation = 0.0;

        for (Eigen::Index j = 0; j < M; ++j) {
            v[j] = K.col(j).transpose() * u;
            vDeviation += pow(aM - vOld[j] * v[j], 2);
        }
        v = aM * v.array().cwiseInverse();

    } while (uDeviation > err && vDeviation > err);

    return K.cwiseProduct(u * v.transpose());
}

double Sinkhorn::distance(const Eigen::MatrixXd& correlationMatrix, double gamma, double eps) {
    return (Pgamma(correlationMatrix, gamma, eps).transpose() * correlationMatrix).trace();
}
