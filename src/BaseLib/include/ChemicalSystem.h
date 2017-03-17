//
// Created by Michael Heuer on 16.03.17.
//

#ifndef AMOLQCGUI_CHEMICALSYSTEM_H
#define AMOLQCGUI_CHEMICALSYSTEM_H

#include <Eigen/Core>

class PositionCollection {
public:

  PositionCollection(const Eigen::VectorXd &pc)
    : positionCollection(pc) {};

  PositionCollection(const unsigned n)
    : positionCollection(Eigen::VectorXd::Zero(n)) {};

  PositionCollection()
    : positionCollection(Eigen::VectorXd::Zero(0)) {};

  Eigen::VectorXd getPositionCollection() {
    return positionCollection;
  };

protected:
  Eigen::VectorXd positionCollection;
};


class NucleusCollection : PositionCollection { // ElementTypeCollection
public:

  //NucleusCollection(PositionCollection pc)
  //  : PositionCollection(pc){}

  NucleusCollection(const Eigen::VectorXd &nc)
    : PositionCollection(nc) {};

  NucleusCollection(const unsigned n)
    : PositionCollection(n) {};

  NucleusCollection()
    : PositionCollection() {};

private:
  unsigned nucleusNumber;
};

class ElectronCollection: PositionCollection  {
  ElectronCollection(const Eigen::VectorXd &ec)
    : PositionCollection(ec) {};

  ElectronCollection(const unsigned n)
    : PositionCollection(n) {};

  ElectronCollection()
    : PositionCollection() {};

private:
  unsigned electronNumber;
};

// can include several molecular geometries or sequences of them
class ChemicalSystem {

};


class MolecularGeometry{

};


#endif //AMOLQCGUI_CHEMICALSYSTEM_H
