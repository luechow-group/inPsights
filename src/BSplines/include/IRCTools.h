//
// Created by Michael Heuer on 08.04.16.
//

#ifndef RTQC_IRCTOOLS_H
#define RTQC_IRCTOOLS_H

namespace Delib {
  class PositionCollection;
  class ElementTypeCollection;
  class MolecularTrajectory;
}

/*! Contains methods to calculate mass-weighted coordinates and distances for Delib types.
 * */
class IRCTools{
public:
  IRCTools();

  Delib::PositionCollection calculateMassWeightedPositionCollection(const Delib::PositionCollection & positions,
                                                                    const Delib::ElementTypeCollection & elements);

  double euclidianDistance(const Delib::PositionCollection &positions1,
                           const Delib::PositionCollection &positions2);

  double calculateMassWeightedDistance(const Delib::PositionCollection & positions1,
                                       const Delib::PositionCollection & positions2,
                                       const Delib::ElementTypeCollection & elements);

  std::vector<double> euclidianDistances(const Delib::MolecularTrajectory & band);

  double length(const Delib::MolecularTrajectory &band);

  std::vector<double> calculateMassWeigtedInterImageDistances(const Delib::MolecularTrajectory & band);

  std::vector<double> calculateDistancesFromIdx(const Delib::MolecularTrajectory & band,
                                                unsigned idx = 0);

  std::vector<double> calculateMassWeightedDistancesFromIdx(const Delib::MolecularTrajectory &band,
                                                            unsigned idx = 0);

  unsigned closestIdx(std::vector<double> const & vec, double value);

  void generateBandWithEqualDistances(Delib::MolecularTrajectory & band,
                                      const unsigned nWantedBeads);

  void generateBandFromSegmentAroundIdx(Delib::MolecularTrajectory & band,
                                                  const unsigned idx,
                                                  const unsigned nNeighbors);

  bool writeMassWeightedDistancesFromIdxToDisk(const std::string &fileName,
                                                   const Delib::MolecularTrajectory &band,
                                                   int idx = -1);

};

#endif //RTQC_IRCTOOLS_H
