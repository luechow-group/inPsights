//
// Created by Michael Heuer on 14.03.18.
//

#ifndef AMOLQCPP_GAUSSKRONRODORDERTYPE_H
#define AMOLQCPP_GAUSSKRONRODORDERTYPE_H

#include "GaussKronrodCartesianIntegration.h"

namespace GaussKronrod{
    using namespace Eigen;
    enum OrderType
    {
        /*A n-point Gauss rule is extende by n+1 Kronrod points given a total of 2n+1 points. */
        GK015 =  22, //   7 point Gauss rule,   8 Kronrod points
        GK021 =  64, //  10 point Gauss rule,  11 Kronrod points
        GK031 =  94, //  15 point Gauss rule,  16 Kronrod points
        GK041 = 124, //  20 point Gauss rule,  20 Kronrod points
        GK051 = 154, //  25 point Gauss rule,  26 Kronrod points
        GK061 = 184, //  30 point Gauss rule,  31 Kronrod points
        GK071 = 214, //  35 point Gauss rule,  36 Kronrod points
        GK081 = 244, //  40 point Gauss rule,  41 Kronrod points
        GK091 = 274, //  45 point Gauss rule,  46 Kronrod points
        GK101 = 304, //  50 point Gauss rule,  51 Kronrod points
        GK121 = 364, //  60 point Gauss rule,  61 Kronrod points
        GK201 = 604, // 100 point Gauss rule, 101 Kronrod points

        NumberOfOrders = 12,
        NotAvailable = 0
    };

    // The n of the nodes of K(2n+1) Gauss-Kronrod rule coincide with those of the
    // n-point Gaussian quadrature rule G(n) for the same measure.
    // K(2n+1)f = If (Integral of f),
    // whenever f is a polynomial of degree less than or equal to 3n+1.

    // Array of all orderTypes to allow iteration
    const std::array<OrderType, static_cast<unsigned>(OrderType::NumberOfOrders)> allOrders = {
            OrderType::GK015, OrderType::GK021, OrderType::GK031, OrderType::GK041, OrderType::GK051,
            OrderType::GK061, OrderType::GK071, OrderType::GK081, OrderType::GK091, OrderType::GK101,
            OrderType::GK121, OrderType::GK201
    };

    static Integrator<double>::QuadratureRule convertOrderTypeToQuadratureRule(const OrderType& orderType){
        switch (orderType){
            case GK015:
                return Integrator<double>::QuadratureRule::GaussKronrod15;
            case GK021:
                return Integrator<double>::QuadratureRule::GaussKronrod21;
            case GK031:
                return Integrator<double>::QuadratureRule::GaussKronrod31;
            case GK041:
                return Integrator<double>::QuadratureRule::GaussKronrod41;
            case GK051:
                return Integrator<double>::QuadratureRule::GaussKronrod51;
            case GK061:
                return Integrator<double>::QuadratureRule::GaussKronrod61;
            case GK071:
                return Integrator<double>::QuadratureRule::GaussKronrod71;
            case GK081:
                return Integrator<double>::QuadratureRule::GaussKronrod81;
            case GK091:
                return Integrator<double>::QuadratureRule::GaussKronrod91;
            case GK101:
                return Integrator<double>::QuadratureRule::GaussKronrod101;
            case GK121:
                return Integrator<double>::QuadratureRule::GaussKronrod121;
            case GK201:
                return Integrator<double>::QuadratureRule::GaussKronrod201;

            default:
                return Integrator<double>::QuadratureRule::GaussKronrod15;
        }
    };

    static Eigen::Integrator<double>::QuadratureRule findAdequateRule(unsigned requiredPolynomialOrder) {
        for (const auto &order : allOrders) {
            if (static_cast<unsigned>(order) >= requiredPolynomialOrder)
                return convertOrderTypeToQuadratureRule(order);
        }
    };
}

#endif //AMOLQCPP_GAUSSKRONRODORDERTYPE_H
