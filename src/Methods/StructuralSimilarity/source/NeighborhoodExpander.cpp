//
// Created by Michael Heuer on 15.05.18.
//

#include "NeighborhoodExpander.h"
#include "CutoffFunction.h"
#include "AngularBasis.h"
#include "Environment.h"
#include <iomanip>
#include <Type.h>

NeighborhoodExpander::NeighborhoodExpander()
        : radialGaussianBasis_(){}

std::complex<double> NeighborhoodExpander::coefficient(double centerToNeighborDistance,
                                                       double theta, double phi,
                                                       double weight, double weightScale,
                                                       double neighborSigma,
                                                       unsigned n, unsigned l, int m) const {
    auto coefficient = radialGaussianBasis_.computeCoefficient(n, l, centerToNeighborDistance, neighborSigma)
                       * AngularBasis::computeCoefficient(l, m, theta, phi);

    return coefficient* weight*weightScale;
}

NeighborhoodExpansion NeighborhoodExpander::expandEnvironment(const Environment& e, int expansionTypeId) const {

    //std::cout << e.atoms_[e.centerId_].toString() << "->\n";//TODO delete

    NeighborhoodExpansion neighborhoodExpansion;

    for (unsigned j = 0; j < unsigned(e.molecularGeometry_.numberOfEntities()); ++j) {
        const auto& center = e.center_;
        const auto& neighbor = e.molecularGeometry_[j];


        //Only add type specific neighbors or ignore the type
        if( int(neighbor.type()) != expansionTypeId && expansionTypeId != int(Type::None)) continue;

        double neighborSigma = 0.5;//TODO arbitrary? change this // double neighborSigma = s_.radial.sigmaAtom;

        Eigen::Vector3d centerToNeighborVector = (neighbor.position()-center);
        double centerToNeighborDistance = centerToNeighborVector.norm();

        // skip this iteration if particle i is outside the cutoff radius
        if (!CutoffFunction::withinCutoffRadiusQ(centerToNeighborDistance)){
            //std::cout << "x";//TODO delete
            //std::cout << "\t" << neighbor.toString() << std::endl;//TODO delete
        } else {
            double weight = 1;
            double weightScale = CutoffFunction::getWeight(centerToNeighborDistance);

            if (centerToNeighborDistance <= ZeroLimits::radiusZero) {//neighbor.position().isApprox(e.center_)) {//TODO FISHY to test for center!
                weight *= ExpansionSettings::Cutoff::centerWeight;
                //std::cout << "*";//TODO delete
            }
            //std::cout << "\t" << neighbor.toString() << "\t";//TODO delete

            double theta, phi;
            if (centerToNeighborDistance > 0.) {
                BoostSphericalHarmonics::ToSphericalCoords(centerToNeighborVector.normalized(), theta, phi);
                if (phi < 0.) phi += 2 * M_PI;
            } else { // center and neighbor positions are identical
                theta = 0.;
                phi = 0.;
                //return something here?
            }

            for (unsigned n = 1; n <= ExpansionSettings::Radial::nmax; ++n) {
                for (unsigned l = 0; l <= ExpansionSettings::Angular::lmax; ++l) {
                    for (int m = -int(l); m < int(l); ++m) {

                        auto coefficient =
                                radialGaussianBasis_.computeCoefficient(n, l, centerToNeighborDistance, neighborSigma)
                                * AngularBasis::computeCoefficient(l, m, theta, phi);
                        //std::cout << AngularBasis::computeCoefficient(l, m, theta, phi) << ", ";
                        coefficient *= weight * weightScale;

                        //std::cout<< std::setprecision(4) << coefficient << ",";
                        neighborhoodExpansion.storeCoefficient(n,l,m,coefficient);
                    }
                }
            }
            //std::cout << std::endl;//TODO delete
        }
    }
    return neighborhoodExpansion;
}
