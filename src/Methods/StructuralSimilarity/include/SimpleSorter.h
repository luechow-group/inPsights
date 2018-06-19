//
// Created by Michael Heuer on 19.06.18.
//

#ifndef AMOLQCPP_SIMPLESORTER_H
#define AMOLQCPP_SIMPLESORTER_H

#include <vector>
#include "MolecularSpectrum.h"
#include "StructuralSimilarity.h"

class SimpleSorter{
public:

    std::vector<std::vector<unsigned >> sort(std::vector<MolecularSpectrum> spectra, double threshold = 0.98){
        assert(!spectra.empty() && "The vector cannot be empty");

        std::vector<std::vector<unsigned >> clusters({0});

        unsigned n = spectra.size();
        for (unsigned i = 1; i < n; ++i) {//iterate over all structures
            std::vector<std::vector<unsigned >>::iterator it;
            for (it = clusters.begin(); it != clusters.end(); ++it){
                unsigned firstElementSpectrumIdx = (*it)[0]; // compare with the spectrum of the first element
                //auto kdist = StructuralSimilarity::kernelDistance(spectra[i],spectra[spectrumIdxOfTheFirstElement]);
                auto kdist = StructuralSimilarity::kernel(spectra[i],spectra[firstElementSpectrumIdx]);
                if(kdist >= threshold) {
                    (*it).push_back(i);
                    break;
                }
            }
            // didn't fit into any other cluster
            if(it == clusters.end()) {
                clusters.push_back({i});
            }
        }

        for (auto& cluster : clusters){
            for(auto& element : cluster){
                std::cout << element << ",";
            }
            std::cout << std::endl;
        }
        return clusters;
    };
};



#endif //AMOLQCPP_SIMPLESORTER_H
