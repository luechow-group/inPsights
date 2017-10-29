//
// Created by heuer on 24.05.17.
//


#include "MolecularGeometry3D.h"
#include "Atom3D.h"
#include "Bond3D.h"

MolecularGeometry3D::MolecularGeometry3D(Qt3DCore::QEntity *root, AtomCollection atomCollection) {

  std::vector<Atom> atoms = atomCollection.atoms();
  std::vector<Atom3D> atoms3D;

  // Draw atoms
  for(std::vector<Atom>::const_iterator atomIt = atoms.begin(); atomIt != atoms.end(); ++ atomIt){
    Eigen::Vector3d vec= (*atomIt).position();
    atoms3D.emplace_back(Atom3D(root,QVector3D(float(vec[0]),float(vec[1]),float(vec[2])), (*atomIt).elementType()));
  }

  float bondDrawingLimit = 1.60*1e-10/AU::length;
  // Draw bonds
  for(std::vector<Atom3D>::const_iterator it1 = atoms3D.begin(); it1 != atoms3D.end(); ++ it1){
    for(std::vector<Atom3D>::const_iterator it2 = it1+1; it2 != atoms3D.end(); ++ it2){
      if( ((*it1).getLocation()-(*it2).getLocation()).length() < bondDrawingLimit) {
        Bond3D(*it1,*it2);
      }
    }
  }
}
