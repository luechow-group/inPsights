//
// Created by Michael Heuer on 10.07.16.
//

#ifndef RTQC_NUMERICALINTEGRATION_H
#define RTQC_NUMERICALINTEGRATION_H

/*! Contains a static function allowing for numberical integration with Simpson's rule*/
class NumericalIntegration{
public:

  //integration based on the Simpson's rule
  template<typename Func>
  static double integrate(Func f, double a, double b, int steps)
  {
    double s = 0;
    double h = (b-a)/steps;
    for (int i = 0; i < steps; ++i) {
      double x = a + h*i;
      s += (f(x) + 4 * f(x + h / 2) + f(x + h)) / 6;
    }
    return h*s;
  }


};

#endif //RTQC_NUMERICALINTEGRATION_H
