//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include <iomanip>
#include <string>

#include "ChemicalSystem.h"
#include "WaveFunctionParser.h"

int main(int argc, char const *argv[]) {

  std::string filename = "t.wf";

  WaveFunctionParser waveFunctionParser(filename);
  waveFunctionParser.readNuclei();

  auto ac = waveFunctionParser.getAtomCollection();

  std::cout << ac.asEigenVector() << std::endl;
}
