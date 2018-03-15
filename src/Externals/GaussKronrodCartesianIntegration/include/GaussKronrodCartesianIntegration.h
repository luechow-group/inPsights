//
// Created by Michael Heuer on 12.03.18.
//

#ifndef AMOLQCPP_GAUSSKRONRODCARTESIANINTEGRATION_H
#define AMOLQCPP_GAUSSKRONRODCARTESIANINTEGRATION_H

/*
 * This file creates a single header to include the NumericalIntegration library.
 * Please note that the order of the includes is important.
 *
 * Warning: MPFR (mpreal.h) is under GPL3 License! This would force us to use GPL3 as well.
 * Thus, we do not want to include it in the numerical integration library.
 */
#include <Eigen/Eigenvalues>
#include <iostream>
//#include "mpreal.h"
#include "NumericalIntegration/Piessens.h"
#include "NumericalIntegration/Monegato.h"
#include "NumericalIntegration/LaurieGautschi.h"
#include "NumericalIntegration/GaussKronrodNodesWeights.h"
#include "NumericalIntegration/ComputeGaussKronrodNodesWeights.h"
#include "NumericalIntegration/Integrator.h"

#include "OrderType.h"

#endif //AMOLQCPP_GAUSSKRONRODCARTESIANINTEGRATION_H
