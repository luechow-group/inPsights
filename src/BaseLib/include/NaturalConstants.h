#ifndef NATURALCONSTANTS_H
#define NATURALCONSTANTS_H

/*! \file DelibConstants.h
 * This header file defines constants commonly used in computational chemistry.
 * Source: http://physics.nist.gov/cuu/Constants/Table/allascii.txt, 03.11.2015.
 */

#include <math.h>

namespace Constant {
    // values from wikipedia // TODO check at NIST gov
    const double pi = M_PI;
    const double boltzmann            = 1.3806488E-23;    // k_B [J/K]
    const double avogadro             = 6.02214129E+23;   // N_A [1/mol]
    const double planckQuantum        = 6.62606957E-34;   // h [J/s]
    const double vacuumPermittivity   = 8.854187817E-12;  // epsilon_0 [A*s/(V*m)]
    const double vacuumPermeability   = 12.566370614E-7;  // mu_0 [N/A^2]
    const double elementaryCharge     = 1.602176565E-19;  // e [C]
    const double electronMass         = 9.10938291E-31;   // m_e [kg]


    // derived

    const double speedOfLight         = sqrt(1./(vacuumPermittivity*vacuumPermeability)); // c0 [m/s]
    //double speedOfLight         = 299792458; // c0 [m/s]
    const double universalGasConstant = boltzmann*avogadro; // R [J/(mol*K)]
    const double reducedPlanckQuantum = planckQuantum/(2*pi); // hbar [J/s]
}

/* Atomic Units */
namespace AU {
    const double length = 4*Constant::pi*Constant::vacuumPermittivity *pow(Constant::reducedPlanckQuantum,2)
                          / (Constant::electronMass*pow(Constant::elementaryCharge,2)); // a0 [m]

    const double energy = Constant::electronMass*pow(Constant::elementaryCharge,4)
                          / (4*Constant::pi*Constant::vacuumPermittivity *pow(Constant::reducedPlanckQuantum,2)); // Eh [J]

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
    const double meter2angstrom     = 1E10;
    const double angstrom2meter     = 1. / meter2angstrom;
    const double meter2bohr         = 1. / AU::length;
    const double bohr2meter         = 1. / meter2bohr;
    const double bohr2angstrom      = AU::length*1E10;
    const double angstrom2bohr      = 1. / bohr2angstrom;

    // conversion between joule, kcal, hartree, electronVolt
    const double joule2hartree      = 1. / AU::energy;
    const double hartree2joule      = 1. / joule2hartree;
    //const double joule2kcal         =
    //const double kcal2joule         =
    const double joule2ev           = 1./ Constant::elementaryCharge;
    const double ev2joule           = 1. / joule2ev;
    const double hartree2ev         = AU::electricPotential;
    const double ev2hartree         = 1. / hartree2ev;
}

#endif // NATURALCONSTANTS_H