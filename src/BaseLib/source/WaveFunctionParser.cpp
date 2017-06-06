#include <sstream>
#include <iostream>
#include <regex>
#include <string>

#include "WaveFunctionParser.h"
#include "ElementInfo.h"

WaveFunctionParser::WaveFunctionParser(const std::string &filename)
        : filename_(filename)
{
  file_.open(filename);
}

WaveFunctionParser::~WaveFunctionParser(){
  file_.close();
}

static std::string extractRegExMatchingSubstring(const std::string str, std::regex rgx){

  std::smatch match;

  if (std::regex_search(str.begin(), str.end(), match, rgx)) {
    std::cout << match[1];
    return match[1];
  }
  else return  "not found";
}


static std::vector<std::string> extractAllRegExMatchingSubstrings(const std::string str, std::regex rgx){
  std::sregex_iterator iter(str.begin(), str.end(), rgx);
  std::sregex_iterator end;

  std::vector<std::string> foundSubstrings;

  while(iter != end)
  {
    for(unsigned i = 0; i < iter->size(); ++i)
    {
      std::string substring = (*iter)[i];
      foundSubstrings.push_back( substring );
    }
    ++iter;
  }

  return foundSubstrings;
}

void WaveFunctionParser::readNuclei() {

  //string str = "File names are readme.txt and my.cmd.";
  //sregex_iterator it(str.begin(), str.end(), reg1);
  //sregex_iterator it_end;

  std::string line;

  std::getline(file_,line);
  std::getline(file_,line);
  std::getline(file_,line);

  std::istringstream iss (line);
  std::string word;

  std::regex chargeNumberRegEx("^charge=(-?[0-9]+),$");
  std::regex spinNumberRegEx("^spin=([0-9]+),?$");
  std::regex angstromRegEx("^geom=angstrom$");
  //std::regex decimalNumberRegEx("^(-?[0-9]+\.[0-9]+)$");

  iss >> word;
  int charge = std::stoi(extractRegExMatchingSubstring(word,chargeNumberRegEx));

  iss >> word;
  auto a =  extractRegExMatchingSubstring(word,spinNumberRegEx);
  std::cout << a;
  unsigned long spin = std::stoul(extractRegExMatchingSubstring(word,spinNumberRegEx));

  iss >> word;
  bool angstromQ = std::regex_match(word,angstromRegEx);

  std::getline(file_,line);
  std::getline(file_,line);
  std::getline(file_,line);
  std::getline(file_,line);
  iss = std::istringstream(line);

  unsigned numberOfNuclei;
  iss >> numberOfNuclei;


  // extract atoms
  Elements::ElementType elementType = Elements::ElementType::none;
  double x;
  double y;
  double z;

  atomCollection_.clear();

  for (unsigned i = 0; i < numberOfNuclei ; ++i) {
    std::getline(file_,line);
    iss = std::istringstream(line);

    iss >> word;
    elementType = Elements::ElementInfo::elementTypeForSymbol(word);

    iss >> x;
    iss >> y;
    iss >> z;

    // if angstrom convert to bohr for amolqc
    if(angstromQ){

      std::cout << AU::length*1e10;
      x /= AU::length*1e10;
      y /= AU::length*1e10;
      z /= AU::length*1e10;
    }

    atomCollection_.addAtom(elementType,x,y,z);
  }
}

