//
// Created by heuer on 19.10.18.
//

#ifndef INPSIGHTS_GLOBALCLUSTERSORTER_H
#define INPSIGHTS_GLOBALCLUSTERSORTER_H

#include "SimilarReferences.h"
#include <DensityBasedScan.h>
#include <Logger.h>

class GlobalClusterSorter {
public:

    GlobalClusterSorter(
            std::vector<Sample> &samples,
            std::vector<SimilarReferences> &globallySimilarMaxima,
            std::vector<std::vector<SimilarReferences>> &globallyClusteredMaxima,
            double similarDistThresh);

    void sort();

private:
    std::vector<Sample> &samples_;
    std::vector<SimilarReferences> &globallySimilarMaxima_;
    std::vector<std::vector<SimilarReferences>> &globallyClusteredMaxima_;
    double similarDistThresh_;
    std::shared_ptr<spdlog::logger> console;


    struct SortElement {
        SortElement(
                std::pair<double, Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch,
                std::vector<SimilarReferences>::iterator it);

        bool operator<(const SortElement &rhs) const;

        std::pair<double, Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch_;
        std::vector<SimilarReferences>::iterator it_;
    };

    static double wrapper(const SimilarReferences &s1, const SimilarReferences &s2);
};

#endif //INPSIGHTS_GLOBALCLUSTERSORTER_H
