//
// Created by Michael Heuer on 10.05.16.
//

#include "BSplineBasis.h"

/*
 * pk is the degree of this basis function of the kth-order derivative spline - not of the zeroth-order spline
 * */
double BSplineBasis::evaluate(const unsigned i,
                              const unsigned pk,
                              const unsigned n,
                              const Eigen::VectorXd &Uk,
                              const double u) const {

  assert( (u >= Uk(pk)) && (u <= Uk(n+1)) && "Parameter has to lie inside of the domain [u_{p},u_{n+1}]");

  double summand1 = 0, summand2 = 0;
  if (pk == 0) {
    /*! the conditional (i == n) &&  (u == Uk(n+1)) closes the last interval [u_{n},u_{n+1}] */
    if ( ((Uk(i) <= u)&&(u < Uk(i+1))) || ((i == n)&&(u == Uk(n+1))) )
      return 1;
    else
      return 0;
  }
  else {
    if (Uk(i+pk) == Uk(i)) summand1 = 0;
    else summand1 = ( u - Uk(i))/(Uk(i+pk) - Uk(i) ) * evaluate(i, pk-1, n, Uk, u);

    if (Uk(i+pk+1) == Uk(i+1)) summand2 = 0;
    else summand2 = ( Uk(i+pk+1) - u) / (Uk(i+pk+1) - Uk(i+1)) * evaluate(i+1, pk-1, n, Uk, u);

    return summand1 + summand2;
  }
}
