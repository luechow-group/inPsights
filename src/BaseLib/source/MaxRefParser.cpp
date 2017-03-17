#include "MaxRefParser.h"
#include "ChemicalSystem.h"

#include <QRegExp>
#include <QFile>
#include <QTextStream>

#include <iostream>


void MaxRefParser::parseFile(const std::string& filename) {
  QFile file(QString::fromStdString(filename));

  if (!file.exists()) return;

  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;

  QTextStream in(&file);

  QString firstline = in.readLine();

  QRegExp regex_firstline("\\s*(\\d+)\\s*");

  int nrAtoms = 0;
  int pos = regex_firstline.indexIn(firstline);
  if (pos > -1) nrAtoms = regex_firstline.cap(1).toInt();

  std::cout << nrAtoms << std::endl;

  QString secondLine = in.readLine();

  int count = 0;
  while (!in.atEnd() && count < nrAtoms) {
    std::string line = in.readLine().toStdString();
    //if (parseLine(line, set)) count++;
  }
}

bool MaxRefParser::parseLine(const std::string& line) {
  QRegExp regex("\\s*([a-zA-Z]{1,3})\\d*\\s+([+-]?\\d*(?:\\.(?:\\d+)?)?)\\s+([+-]?\\d*(?:\\.(?:\\d+)?)?)\\s+([+-]?\\d*(?:\\.(?:\\d+)?)?)\\s*");

  int pos = regex.indexIn(QString::fromStdString(line));

  if (pos > -1) {
    std::string symbol(regex.cap(1).toUtf8());
    //Atom* atom = new Atom(ElementInfo::instance()[symbol].first);
    //auto positionInAngstrom = Eigen::Vector3d(regex.cap(2).toDouble(), regex.cap(3).toDouble(), regex.cap(4)
    //  .toDouble());
    //atom->setPosition(positionInAngstrom * angstrom2bohr);
    //set.addChild(atom);
    return true;
  }
  else
    return false;
}
