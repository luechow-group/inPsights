// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "Sinkhorn.h"

using namespace testing;
using namespace SOAP;

class ASinkhornTest : public ::testing::Test {
public:

    double machineEpsilon = std::numeric_limits<double>::epsilon(); // same as default values
    double comparisonEpsilon = machineEpsilon * 1e2;

    static Eigen::MatrixXd sinkhornSoapxxReference(const Eigen::MatrixXd &matrix,
                                                   double gamma = General::settings.sinkhornGamma(),
                                                   double eps = std::numeric_limits<double>::epsilon()) {
        long nx = matrix.rows();
        long ny = matrix.cols();
        std::vector<double> u(nx), ou(nx), v(ny);
        double ax = 1.0/double(nx), ay=double(1.0)/ny;
        Eigen::MatrixXd Kg(nx,ny);
        for (int i=0; i<nx; ++i) u[i]=1.0;
        for (int i=0; i<ny; ++i) v[i]=1.0;
        double lambda=1.0/gamma, terr=eps*eps, derr;

        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                Kg(i, j) = std::exp(-(1 - matrix(i, j)) * lambda);
        do {
            // u<-1.0/Kg.v
            for (int i = 0; i < nx; ++i) {
                ou[i] = u[i];
                u[i] = 0.0;
            }
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    u[i] += Kg(i, j) * v[j];
            // at this point we can compute how far off unity we are
            derr = 0.0;
            for (int i = 0; i < nx; ++i)
                derr += (ax - ou[i] * u[i]) * (ax - ou[i] * u[i]);
            for (int i = 0; i < nx; ++i)
                u[i] = ax / u[i];

            // v<-1.0/Kg.u
            for (int i = 0; i < ny; ++i)
                v[i] = 0.0;
            for (int i = 0; i < ny; ++i)
                for (int j = 0; j < nx; ++j) v[i] += Kg(j, i) * u[j];
            for (int i = 0; i < ny; ++i)
                v[i] = ay / v[i];
            //std::cerr<<derr<<"\n";

        } while (derr > terr);

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                Kg(i, j) *= u[i] * v[j];
            }
        }

        return Kg;
    }

    static Eigen::MatrixXd sinkhornSoapPythonReference(const Eigen::MatrixXd &matrix,
                                                       double gamma = General::settings.sinkhornGamma()) {
/* Sinkhorn algorithm */
        auto MAX_TOTAL = matrix.rows();
        Eigen::VectorXd v = Eigen::VectorXd::Zero(MAX_TOTAL);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(MAX_TOTAL);
        Eigen::VectorXd en = Eigen::VectorXd::Zero(MAX_TOTAL);
        Eigen::VectorXd temp(MAX_TOTAL);
        Eigen::MatrixXd Sink = Eigen::MatrixXd::Zero(MAX_TOTAL, MAX_TOTAL);
        int ITERATIONS_NO = 20;


        for (int i = 0; i < MAX_TOTAL; i++) {
            for (int j = 0; j < MAX_TOTAL; j++) {
                Sink(i, j) = exp((matrix(i, j) - 1.0) / gamma);
            }
        }

        for (int i = 0; i < MAX_TOTAL; i++) en(i) = 1.0 / double(MAX_TOTAL);
        v = en;
        for (int counter = 0; counter < ITERATIONS_NO; counter++) {
            temp = Sink * v;
            for (int i = 0; i < MAX_TOTAL; i++) u(i) = en(i) / temp(i);
            temp = Sink.transpose() * u;
            for (int i = 0; i < MAX_TOTAL; i++) v(i) = en(i) / temp(i);
        }

        for (int i = 0; i < MAX_TOTAL; i++) {
            for (int j = 0; j < MAX_TOTAL; j++) {
                Sink(i, j) *= u(i) * v(j);
            }
        }
        return Sink;
    }

    static void compareWithReferenceImplementations(const Eigen::MatrixXd &correlationMatrix,
                                             double intendedError = std::numeric_limits<double>::epsilon(),
                                             double comparisonEpsilion = std::numeric_limits<double>::epsilon()) {
        auto gamma = 1.0;
        auto ownImplementation = Sinkhorn::Pgamma(correlationMatrix, gamma, intendedError);
        auto referenceImplementation1 = sinkhornSoapxxReference(correlationMatrix, gamma, intendedError);
        ASSERT_TRUE(ownImplementation.isApprox(referenceImplementation1, comparisonEpsilion));

        if (correlationMatrix.rows() == correlationMatrix.cols()) {
            auto referenceImplementation2 = sinkhornSoapPythonReference(correlationMatrix, gamma);
            ASSERT_TRUE(referenceImplementation1.isApprox(referenceImplementation2, comparisonEpsilion));
        }
    }
};


void doublyStochasticCheck(const Eigen::MatrixXd &regularizedMatrix,
                           double eps = std::numeric_limits<double>::epsilon()) {
    long N = regularizedMatrix.rows();
    long M = regularizedMatrix.cols();

    for (unsigned i = 0; i < N; ++i) {
        ASSERT_NEAR(regularizedMatrix.row(i).sum(), 1. / double(N), eps);
    }

    for (unsigned i = 0; i < M; ++i) {
        ASSERT_NEAR(regularizedMatrix.col(i).sum(), 1. / double(M), eps);
    }

    ASSERT_NEAR(regularizedMatrix.sum(), 1.0, eps);
}

TEST_F(ASinkhornTest, BestMatchOfPermutedIdenticalEnvironments) {
    Eigen::MatrixXd C(2, 2); // covariance matrix of two permuted particles
    C << 0, 1, 1, 0;  // the system is composed of two permuted but otherwise identical environments
    Eigen::MatrixXd Pref(2, 2);
    Pref << 0, 0.5, 0.5, 0;
    Eigen::MatrixXd expected(2, 2);
    expected << 0.5, 0, 0, 0.5;

    compareWithReferenceImplementations(C, machineEpsilon, comparisonEpsilon);
    Eigen::MatrixXd P = Sinkhorn::Pgamma(C, 0.00001, machineEpsilon);
    doublyStochasticCheck(P, comparisonEpsilon);

    Eigen::MatrixXd result = P * C.transpose();

    ASSERT_TRUE(result.isApprox(expected));
    ASSERT_EQ(result.trace(), 1); // the similarity must equal 1
}

TEST_F(ASinkhornTest, DoublyStochasticCheck_2x2) {
    Eigen::MatrixXd C(2, 2);
    C << 0.6, 0.5, 0.5, 0.2;
    compareWithReferenceImplementations(C, machineEpsilon, comparisonEpsilon);
    doublyStochasticCheck(Sinkhorn::Pgamma(C, 1, machineEpsilon), comparisonEpsilon);
}

TEST_F(ASinkhornTest, DoublyStochasticCheck_20x20RandomMatrix) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(20, 20);
    compareWithReferenceImplementations(C, machineEpsilon, comparisonEpsilon);
    doublyStochasticCheck(Sinkhorn::Pgamma(C, 1, machineEpsilon), comparisonEpsilon);
}

TEST_F(ASinkhornTest, DoublyStochasticCheck_20x10RandomMatrix) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(20, 10);
    compareWithReferenceImplementations(C, machineEpsilon, comparisonEpsilon);
    doublyStochasticCheck(Sinkhorn::Pgamma(C, 1, machineEpsilon), comparisonEpsilon);
}

TEST_F(ASinkhornTest, DoublyStochasticCheck_10x20RandomMatrix) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(10, 20);
    compareWithReferenceImplementations(C, machineEpsilon, comparisonEpsilon);
    doublyStochasticCheck(Sinkhorn::Pgamma(C, 1, machineEpsilon), comparisonEpsilon);
}

TEST_F(ASinkhornTest, DoublyStochasticCheck_100x100RandomMatrix) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(100, 100);
    compareWithReferenceImplementations(C, machineEpsilon, comparisonEpsilon);
    doublyStochasticCheck(Sinkhorn::Pgamma(C, 1, machineEpsilon), comparisonEpsilon);
}
