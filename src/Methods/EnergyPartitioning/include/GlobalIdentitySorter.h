//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_GLOBALIDENTITYSORTER_H
#define INPSIGHTS_GLOBALIDENTITYSORTER_H

#include "Reference.h"
#include "Sample.h"
#include <Logger.h>
#include <utility>
#include <vector>

class GlobalIdentiySorter {
public:

    GlobalIdentiySorter(std::vector<Reference> &references, std::vector<Sample> &samples,
                        double distThresh, double increment);
    bool sort();

private:
    void subLoop(
            std::vector<Reference>::iterator &beginIt,
            std::vector<Reference>::iterator &it,
            std::vector<Reference>::iterator &endIt);

    // TODO This method should be located inside of a reference container class
    void addReference(
            const std::vector<Reference>::iterator &beginIt,
            std::vector<Reference>::iterator &it,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const;

    std::vector<Reference> &references_;
    std::vector<Sample> &samples_;
    double increment_, distThresh_;
};


#endif //INPSIGHTS_GLOBALIDENTITYSORTER_H
