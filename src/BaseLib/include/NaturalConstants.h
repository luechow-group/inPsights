#ifndef NATURALCONSTANTS_H
#define NATURALCONSTANTS_H

/* This header file defines constants commonly used in computational chemistry.
 * Source: http://physics.nist.gov/cuu/Constants/Table/allascii.txt, 03.11.2015.
 */

#include <cmath>

namespace Constant {
    // values from https://physics.nist.gov, accessed on march 01, 2018
    const double pi = M_PI;
    const double boltzmann = 1.38064852E-23; // k_B [J/K] (79)
    const double avogadro = 6.022140857E+23; // N_A [1/mol] (74)
    const double planckQuantum = 6.626070040E-34; // h [J/s] (81)
    const double elementaryCharge = 1.6021766208E-19; // e [C] (98)
    const double electronMass = 9.10938356E-31; // m_e [kg] (11)
    const double speedOfLight = 299792458; // c0 [m/s] (exact)
    const double electricConstant = 8.854187817E-12; // epsilon_0 [A*s/(V*m)] (exact)
    const double magneticConstant = 12.566370614E-7; // mu_0 [N/A^2] (exact)
    //const double speedOfLight = std::sqrt(1./(electricConstant*magneticConstant)); // c0 [m/s]

    // derived
    const double universalGasConstant = boltzmann*avogadro; // R [J/(mol*K)]
    const double reducedPlanckQuantum = planckQuantum/(2*pi); // hbar [J/s]
}

/* Atomic Units */
namespace AU {
    const double length = 0.52917721067E-10; // a0 [m] (12)
    //const double length = 4*Constant::pi*Constant::electricConstant *pow(Constant::reducedPlanckQuantum,2)
    //                      / (Constant::electronMass*pow(Constant::elementaryCharge,2)); // a0 [m]
    const double energy = 4.359744650E-18; // Eh [J] (54)
    //const double energy = Constant::electronMass*pow(Constant::elementaryCharge,4)
    //                      / (4*Constant::pi*Constant::electricConstant *pow(Constant::reducedPlanckQuantum,2)); // Eh [J]

    const double time = Constant::reducedPlanckQuantum/energy; // [s]

    const double velocity = length*energy / Constant::reducedPlanckQuantum; // [m/s]

    const double force = energy/length; // [J/m]

    const double temperature = energy/Constant::boltzmann; // [K]

    const double pressure = force/pow(length,2); // [J/m^2]

    const double electricField = energy/(Constant::elementaryCharge*length);// [J/(C*m)]

    const double electricPotential = energy/Constant::elementaryCharge; // [J/C]

    const double dipoleMoment = Constant::elementaryCharge * length; // [C*m]
}

/* Conversion Factors */
namespace ConversionFactors {

    // conversion between meter, angstrom, bohr
    const double meter2angstrom = 1E10;
    const double angstrom2meter = 1. / meter2angstrom;
    const double meter2bohr = 1. / AU::length;
    const double bohr2meter = 1. / meter2bohr;
    const double bohr2angstrom = AU::length*1E10;
    const double angstrom2bohr = 1. / bohr2angstrom;

    // conversion between joule, kcal, hartree, electronVolt
    const double joule2hartree = 1. / AU::energy;
    const double hartree2joule = 1. / joule2hartree;
    //const double joule2kcal = ?
    //const double kcal2joule = ?
    const double joule2ev = 1./ Constant::elementaryCharge;
    const double ev2joule = 1. / joule2ev;
    const double hartree2ev = AU::electricPotential;
    const double ev2hartree = 1. / hartree2ev;

    const double deg2rad = 2.*M_PI/360.;
    const double rad2deg = 1./deg2rad;
}

#endif // NATURALCONSTANTS_H