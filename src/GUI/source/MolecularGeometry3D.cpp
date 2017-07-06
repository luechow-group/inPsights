//
// Created by heuer on 24.05.17.
//

#include "MolecularGeometry3D.h"
#include "Atom3D.h"

MolecularGeometry3D::MolecularGeometry3D(Qt3DCore::QEntity *root, AtomCollection atomCollection) {

  std::vector<Atom> atoms = atomCollection.atoms();

  for(std::vector<Atom>::const_iterator atomIt = atoms.begin(); atomIt != atoms.end(); ++ atomIt){
    Eigen::Vector3d vec= (*atomIt).coordinates();
    Atom3D(root,QVector3D(float(vec[0]),float(vec[1]),float(vec[2])), (*atomIt).elementType());
  }
}
