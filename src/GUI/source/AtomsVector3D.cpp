//
// Created by heuer on 24.05.17.
//


#include "AtomsVector3D.h"
#include "Atom3D.h"
#include "Bond3D.h"

AtomsVector3D::AtomsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector) {

  /*TODO Refactor*/
  std::vector<Atom3D> atoms3D;

  // Draw atoms
  for (long i = 0; i < atomsVector.numberOfEntities(); ++i) {
    Eigen::Vector3d vec= atomsVector[i].position();
    atoms3D.emplace_back(Atom3D(root,
                                QVector3D(float(vec[0]),float(vec[1]),float(vec[2])),
                                atomsVector.typesVector()[i]));
  }

  // Draw bonds
  auto bondDrawingLimit = float(1.50*1e-10/AU::length);
  for(std::vector<Atom3D>::const_iterator it1 = atoms3D.begin(); it1 != atoms3D.end(); ++ it1){
    for(auto it2 = it1+1; it2 != atoms3D.end(); ++ it2){
      if( ((*it1).getLocation()-(*it2).getLocation()).length() < bondDrawingLimit) {
        Bond3D(*it1,*it2);
      }
    }
  }
}
