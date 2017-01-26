//
// Created by Michael Heuer on 26.01.17.
//

#ifndef AMOLQCGUI_NUMERCIALDIFFERENTIATOR_H
#define AMOLQCGUI_NUMERCIALDIFFERENTIATOR_H

namespace Eigen {
  class VectorXd;
};

class NumericalDifferentiator{
public:
  NumericalDifferentiator();

  Eigen::VectorXd differentiate(double (*function)(Eigen::VectorXd));

  double partialDerivative(double (*function)(Eigen::VectorXd), unsigned dim);


  unsigned calculateOptimalDifference(double a, double b);

  double centralDifferences(double (*f)(double),double x0,double h,int N,int diff_degree);

private:

};

#endif //AMOLQCGUI_NUMERCIALDIFFERENTIATOR_H
