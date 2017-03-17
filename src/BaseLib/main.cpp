//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include <iomanip>
#include <string>

#include "ChemicalSystem.h"
#include "MaxRefParser.h"

int main(int argc, char const *argv[]) {

  Eigen::VectorXd t(3);
  t << 1,2,3;

  PositionCollection pc;

  std::string filename;
  filename = "./Ethane-max.ref";

  MaxRefParser maxRefParser;
  maxRefParser.parseFile(filename);


  std::cout << pc.getPositionCollection().size() << std::endl;
  std::cout << pc.getPositionCollection() << std::endl;

}
