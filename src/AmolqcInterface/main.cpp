#include <iostream>
#include <iomanip>
#include "ElectronicWaveFunction.h"

int main(int argc, char const *argv[]) {

  //TODO Make one global resources folder.
  auto &wf = ElectronicWaveFunction::getInstance("Ethane-em-5.wf");

  std::cout << wf.getAtomCollection() << std::endl;

  return 0;
}
