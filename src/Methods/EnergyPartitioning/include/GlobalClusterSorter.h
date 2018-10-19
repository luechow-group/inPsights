//
// Created by heuer on 19.10.18.
//

#ifndef AMOLQCPP_GLOBALCLUSTERSORTER_H
#define AMOLQCPP_GLOBALCLUSTERSORTER_H

#include "Reference.h"
#include <DensityBasedScan.h>
#include <Logger.h>


double wrapper(const SimilarReferences& s1, const SimilarReferences& s2) {
    //TODO add value to the metric
    return Metrics::bestMatchNorm<Eigen::Infinity,2>(
            s1.representativeReference().maximum(),
            s2.representativeReference().maximum());
};

struct SortElement{
    SortElement(
            std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch,
            std::vector<SimilarReferences>::iterator it)
            : bestMatch_(std::move(bestMatch)), it_(it)
    {}

    bool operator < (const SortElement& rhs) const {
        return bestMatch_.first < rhs.bestMatch_.first;
    }

    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch_;
    std::vector<SimilarReferences>::iterator it_;
};

class GlobalClusterSorter{
public:

    GlobalClusterSorter(
            std::vector<SimilarReferences>& globallySimilarMaxima,
            std::vector<std::vector<SimilarReferences>>& globallyClusteredMaxima,
            double similarDistThresh )
    :
    globallySimilarMaxima_(globallySimilarMaxima),
    globallyClusteredMaxima_(globallyClusteredMaxima),
    similarDistThresh_(similarDistThresh),
    console(spdlog::get(Logger::name))
    {
        if(!console){
            Logger::initialize();
            console = spdlog::get(Logger::name);
        };
    }

    void sort(){
        DensityBasedScan<double, SimilarReferences, wrapper> dbscan(globallySimilarMaxima_);
        auto nClusters = dbscan.findClusters(similarDistThresh_*2+0.01, 1); // why multiplication by 2 is needed?
        //TODO Permutations must be stored! Own DBSCAN implementation?
        auto labels = dbscan.getLabels();
        console->info("number of clusters {}",nClusters);

        globallyClusteredMaxima_.resize(nClusters);

        for (int i = 0; i < nClusters; ++i) {
            for (auto it = globallySimilarMaxima_.begin(); it != globallySimilarMaxima_.end(); ++it) {
                auto label = labels[std::distance(globallySimilarMaxima_.begin(),it)];
                if (label == i) {
                    globallyClusteredMaxima_[i].emplace_back(std::move(*it));
                }
            }
        }

        // order clusters by best match distance
        for (auto& cluster : globallyClusteredMaxima_) {
            std::sort(cluster.begin(), cluster.end());

            for (auto i = cluster.begin(); i != cluster.end(); ++i) {
                std::vector<SortElement> bestMatchDistances;

                for (auto j = i + 1; j != cluster.end(); ++j) {
                    bestMatchDistances.emplace_back(
                            SortElement(Metrics::bestMatch<Eigen::Infinity, 2>(
                                    j.base()->representativeReference().maximum(),
                                    i.base()->representativeReference().maximum()), j)
                    );
                }

                if (bestMatchDistances.size() > 1) {
                    // find the SimilarReferences object whose representativeReference is closest to the i
                    auto minIt = std::min_element(bestMatchDistances.begin(), bestMatchDistances.end());
                    // permute and swap
                    if (i + 1 != minIt.base()->it_)
                        std::iter_swap(i + 1, minIt.base()->it_);
                }
            }
        }
    }


private:
    std::vector<SimilarReferences>& globallySimilarMaxima_;
    std::vector<std::vector<SimilarReferences>>& globallyClusteredMaxima_;
    double similarDistThresh_;
    std::shared_ptr<spdlog::logger> console;
};

#endif //AMOLQCPP_GLOBALCLUSTERSORTER_H
