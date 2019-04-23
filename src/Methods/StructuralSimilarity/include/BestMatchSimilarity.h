//
// Created by heuer on 03.04.19.
//

#ifndef INPSIGHTS_BESTMATCHSIMILARITY_H
#define INPSIGHTS_BESTMATCHSIMILARITY_H

#include "BestMatch.h"
#include <MolecularSpectrum.h>

namespace BestMatch {
    namespace Similarity {
        Result compare(
                MolecularGeometry permutee, const MolecularGeometry &reference,
                bool spinSpecificQ = false, bool flipSpinsQ = false);

        Result compare(const SOAP::MolecularSpectrum &permutee, const SOAP::MolecularSpectrum &reference);
    }
}

#endif //INPSIGHTS_BESTMATCHSIMILARITY_H
