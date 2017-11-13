#include <sstream>
#include <iostream>
#include <regex>
#include <string>
#include <algorithm>

#include "WfFileImporter.h"
#include "ElementInfo.h"

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
      else if ( splittedElement[0] == "geom") angstromQ_ = splittedElement[1] == "angstrom";
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

AtomCollection WfFileImporter::getAtomCollection() {

  AtomCollection atomCollection;
  auto startLine = findTag("$geom", 0).second+1;
  numberOfNuclei_ = std::stoul(strip(getLine(startLine)));

    for (unsigned i = 1; i <= numberOfNuclei_; ++i) {
      std::vector<std::string> lineElements = split(getLine(startLine+i));
      Elements::ElementType elementType = Elements::ElementInfo::elementTypeForSymbol(lineElements[0]);
      double x = std::stod(lineElements[1]);
      double y = std::stod(lineElements[2]);
      double z = std::stod(lineElements[3]);

      // if angstrom convert to bohr for amolqc
      if(angstromQ_){
        x /= AU::length*1e10;
        y /= AU::length*1e10;
        z /= AU::length*1e10;
      }
      atomCollection.addAtom(x,y,z,elementType);
    }
    return atomCollection;
}
