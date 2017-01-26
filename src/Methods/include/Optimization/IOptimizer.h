//
// Created by Michael Heuer on 23.01.17.
//

#ifndef AMOLQCGUI_IOPTIMIZER_H
#define AMOLQCGUI_IOPTIMIZER_H

#include <vector>
#include <Eigen/Core>

/*struct CartesionCoordinate{
  double x,y,z;
};*/

class IOptimizer {
public:
  IOptimizer();


  virtual ~IOptimizer() {};
  virtual void performStep() = 0;

  void fetchGradient();
  void iterate();
  bool converged();


protected:
  unsigned stepCount_;
  Eigen::VectorXd positions_, gradients_;

  //std::vector<CartesionCoordinate> positions_, gradients_;
};

#endif //AMOLQCGUI_IOPTIMIZER_H
