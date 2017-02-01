//
// Created by Michael Heuer on 26.01.17.
//

#include "NumercialDifferentiator.h"
#include <unsupported/Eigen/NumericalDiff>

#include "iostream"

NumericalDifferentiator::NumericalDifferentiator() {

}


Eigen::VectorXd NumericalDifferentiator::differentiate(double (*function)(Eigen::VectorXd)) {
  return Eigen::VectorXd();
}


double NumericalDifferentiator::partialDerivative(double (*function)(Eigen::VectorXd), unsigned dim) {

  Eigen::VectorXd vec;
  double v = function(vec);



  //partialDerivative();
  return 0;
}

unsigned NumericalDifferentiator::calculateOptimalDifference(double a, double b) {
  return 0;
}

double NumericalDifferentiator::centralDifferences(double (*f)(double), double x0, double h, int N, int diff_degree) {
  if (diff_degree == 1) {
    if (N == 3) {
      return (f(x0 + h) - f(x0 - h)) / (2 * h);
    } else if (N == 5) {
      return (f(x0 - 2 * h) - 8 * f(x0 - h) + 8 * f(x0 + h) - f(x0 + 2 * h)) / (12 * h);
    } else if (N == 7) {
      return
        (-f(x0 - 3 * h) + 9 * f(x0 - 2 * h) - 45 * f(x0 - h) + 45 * f(x0 + h) - 9 * f(x0 + 2 * h) + f(x0 + 3 * h)) /
        (60 * h);
    } else if (N == 9) {
      return (3 * f(x0 - 4 * h) - 32 * f(x0 - 3 * h) + 168 * f(x0 - 2 * h) - 672 * f(x0 - h) + 672 * f(x0 + h)
              - 168 * f(x0 + 2 * h) + 32 * f(x0 + 3 * h) - 3 * f(x0 + 4 * h)) / (840 * h);
    } else {
      //logs("s","error - no differentiation formula implemented for this number of points!");
      exit(0);
    }
  } else {
    if (N == 3) {
      return (centralDifferences(f, x0 + h, h, N, diff_degree - 1)
              - centralDifferences(f, x0 - h, h, N, diff_degree - 1)) / (2 * h);
    } else if (N == 5) {
      return (centralDifferences(f, x0 - 2 * h, h, N, diff_degree - 1)
              - 8 * centralDifferences(f, x0 - h, h, N, diff_degree - 1)
              + 8 * centralDifferences(f, x0 + h, h, N, diff_degree - 1)
              - centralDifferences(f, x0 + 2 * h, h, N, diff_degree - 1)) / (12 * h);
    } else if (N == 7) {
      return (-centralDifferences(f, x0 - 3 * h, h, N, diff_degree - 1)
              + 9 * centralDifferences(f, x0 - 2 * h, h, N, diff_degree - 1)
              - 45 * centralDifferences(f, x0 - h, h, N, diff_degree - 1)
              + 45 * centralDifferences(f, x0 + h, h, N, diff_degree - 1)
              - 9 * centralDifferences(f, x0 + 2 * h, h, N, diff_degree - 1)
              + centralDifferences(f, x0 + 3 * h, h, N, diff_degree - 1)) / (60 * h);
    } else if (N == 9) {
      return (3 * centralDifferences(f, x0 - 4 * h, h, N, diff_degree - 1)
              - 32 * centralDifferences(f, x0 - 3 * h, h, N, diff_degree - 1)
              + 168 * centralDifferences(f, x0 - 2 * h, h, N, diff_degree - 1)
              - 672 * centralDifferences(f, x0 - h, h, N, diff_degree - 1)
              + 672 * centralDifferences(f, x0 + h, h, N, diff_degree - 1)
              - 168 * centralDifferences(f, x0 + 2 * h, h, N, diff_degree - 1)
              + 32 * centralDifferences(f, x0 + 3 * h, h, N, diff_degree - 1)
              - 3 * centralDifferences(f, x0 + 4 * h, h, N, diff_degree - 1)) / (840 * h);
    } else {
      //logs("s","error - no differentiation formula implemented for this number of points!");
      exit(0);
    }
  }
}
