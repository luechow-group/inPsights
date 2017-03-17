//
// Created by Michael Heuer on 09.01.17.
//

#ifndef MAXREFPARSER_H
#define MAXREFPARSER_H

#include <string>

/*! Class for the conversion of ".xyz" files to AtomSet.
 * The parsing relies on Qt and the conversion to doubles works also when numbers in the current locale are of the form "x,xxx".
 */
class MaxRefParser {
public:
  static void parseFile(const std::string& filename); //, AtomSet& set);

  static bool parseLine(const std::string& line);//, AtomSet& set);
};

#endif //MAXREFPARSER_H
