//
// Created by heuer on 24.05.17.
//

#ifndef AMOLQCGUI_ATOMCOLLECTION_H
#define AMOLQCGUI_ATOMCOLLECTION_H

#include "Atom.h"

class AtomCollection{
public:
    AtomCollection(){};

    void addAtom(const Elements::ElementType& elementType,
                 const double x, const double y, const double z);
    
    void addAtom(const Elements::ElementType& elementType,
                 const Eigen::Vector3d& coordinates );
    Atom atom(const unsigned index){
      return atoms_[index];
    };

    std::vector<Atom> atoms(){ return atoms_; };


    void clear(){
      atoms_.clear();
    }

    Eigen::VectorXd asEigenVector();

private:
    std::vector<Atom> atoms_;
};

#endif //AMOLQCGUI_ATOMCOLLECTION_H
