//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_ENVIRONMENTEXPANDER_H
#define AMOLQCPP_ENVIRONMENTEXPANDER_H

#include "RadialGaussianBasis.h"
#include "NeighborhoodExpansion.h"
#include "ElementType.h"

class Environment;

class NeighborhoodExpander{
public:
    explicit NeighborhoodExpander();

    std::complex<double> coefficient(double centerToNeighborDistance,
                                     double theta, double phi,
                                     double weight, double weightScale,
                                     double neighborSigma,
                                     unsigned n, unsigned l, int m) const;

    NeighborhoodExpansion expandEnvironment(const Environment& e,
                                            Elements::ElementType expansionType = Elements::ElementType::none) const;

private:
    RadialGaussianBasis radialGaussianBasis_;
};

#endif //AMOLQCPP_ENVIRONMENTEXPANDER_H
