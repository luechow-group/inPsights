
#include "AtomCollection.h"

void AtomCollection::addAtom(const Elements::ElementType &elementType,
                             const double x, const double y, const double z) {
  atoms_.emplace_back(Atom(elementType, Eigen::Vector3d(x,y,z)));
};

void AtomCollection::addAtom(const Elements::ElementType &elementType,
                             const Eigen::Vector3d &coordinates) {
  atoms_.emplace_back(Atom(elementType, coordinates));
};

Eigen::VectorXd AtomCollection::asEigenVector() {

  Eigen::VectorXd coordinates(atoms_.size() * 3);

  for (std::vector<Atom>::const_iterator it = atoms_.begin(); it != atoms_.end(); ++it) {
    coordinates.segment((it - atoms_.begin())*3, 3) = (*it).coordinates();
  }

  return coordinates;
}
