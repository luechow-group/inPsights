//
// Created by Michael Heuer on 18.05.18.
//

#include "Sinkhorn.h"

Eigen::MatrixXd Sinkhorn::Pgamma(const Eigen::MatrixXd &C, double gamma, double eps){
    auto N = unsigned(C.rows());
    auto M = unsigned(C.cols());


    Eigen::VectorXd u = Eigen::VectorXd::Constant(N,1./double(N));
    Eigen::VectorXd v = Eigen::VectorXd::Constant(M,1./double(M));
    Eigen::VectorXd uOld(N),vOld(M);
    Eigen::MatrixXd K(N,M);

    double aN = 1.0/double(N), aM=1.0/double(M);
    double uDeviation, vDeviation;
    double err = pow(eps,2);

    double lambda = 1./gamma;
    K = (lambda*(C-Eigen::MatrixXd::Ones(N,M))).array().exp();

    do {
        uOld = u;
        uDeviation = 0.0;

        for (unsigned i=0; i<N; ++i) {
            u[i] = K.row(i) * v;
            // at this point we can compute how far off unity we are
            uDeviation += pow(aN - uOld[i]*u[i],2);
        }
        u = aN*u.array().cwiseInverse();

        vOld = v;
        vDeviation = 0.0;

        for (unsigned j=0; j<M; ++j) {
            v[j] = K.col(j).transpose() * u;
            vDeviation += pow(aM -  vOld[j]*v[j],2);
        }
        v = aM*v.array().cwiseInverse();

    } while ( uDeviation > err && vDeviation > err);

    return K.cwiseProduct(u*v.transpose());

}

double Sinkhorn::distance(Eigen::MatrixXd correlationMatrix, double gamma, double eps) {
    return (Pgamma(correlationMatrix,gamma,eps).transpose()*correlationMatrix).trace();
}
