
#include <WfFileImporter.h>
#include <ElementInfo.h>
#include <TypesVector.h>

WfFileImporter::WfFileImporter(const std::string &filename)
        : AmolqcImporter(filename)
{
  readGeneralBlock();
}

void WfFileImporter::readGeneralBlock() {
  auto startLine = findTag("$general", 0).second;
  auto endLine = findTag("$end", startLine).second;

  for (unsigned long i = startLine+1; i < endLine; ++i) {
    auto lineElements = split(strip(getLine(i)), ',');
    for (const auto &elem : lineElements){

      auto splittedElement = split(elem, '=');

      if ( splittedElement[0] == "title") title_ = splittedElement[1];
      else if ( splittedElement[0] == "basis") basis_ = splittedElement[1];
      else if ( splittedElement[0] == "jastrow") jastrow_ = splittedElement[1];
      else if ( splittedElement[0] == "geom") bohrQ_ = splittedElement[1] == "bohr"; // TODO Implement better detection
      else if ( splittedElement[0] == "charge") charge_ = std::stoi(splittedElement[1]);
      else if ( splittedElement[0] == "spin") multiplicity_ = std::stoul(splittedElement[1]);
    }
  }
}

std::pair<bool, unsigned long> WfFileImporter::findTag(const std::string &tag, unsigned long startLine) {
  unsigned long lineIdx = startLine;

  while( (strip(getLine(lineIdx)).empty() || split(getLine(lineIdx))[0] != tag)
         && lineIdx < lines_.size()){
    lineIdx++;
  }
  return std::make_pair(lineIdx < lines_.size(), lineIdx);
}

AtomsVector WfFileImporter::getAtomsVector() {

  AtomsVector atomsVector;
  auto startLine = findTag("$geom", 0).second+1;
  numberOfNuclei_ = std::stoul(strip(getLine(startLine)));

    for (unsigned i = 1; i <= numberOfNuclei_; ++i) {
      std::vector<std::string> lineElements = split(getLine(startLine+i));
      Element elementType = Elements::ElementInfo::elementTypeFromSymbol(lineElements[0]);
      double x = std::stod(lineElements[1]);
      double y = std::stod(lineElements[2]);
      double z = std::stod(lineElements[3]);

      // The internal amolqc and inPsights default is bohr, if not bohr, convert to bohr
      if(!bohrQ_){
        x *= ConversionFactors::angstrom2bohr;
        y *= ConversionFactors::angstrom2bohr;
        z *= ConversionFactors::angstrom2bohr;
      }
      atomsVector.append({elementType,{x,y,z}});
    }
  return atomsVector;
}

unsigned long WfFileImporter::getNumberOfElectrons() {

  ElementTypesVector etc = getAtomsVector().typesVector();
  unsigned long numberOfElectrons = 0;

  for (unsigned long i = 0; i < etc.numberOfEntities(); ++i) {
    numberOfElectrons += Elements::ElementInfo::Z(etc[i]);
  }
  numberOfElectrons -= charge_;

  return numberOfElectrons;
}

unsigned long WfFileImporter::getNumberOfAlphaElectrons() {
  return (getNumberOfElectrons()+(getMultiplicity()-1))/2;
}

unsigned long WfFileImporter::getNumberOfBetaElectrons() {
  return getNumberOfElectrons()-getNumberOfAlphaElectrons();
}
