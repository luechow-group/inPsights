//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "Sample.h"
#include "Logger.h"

#include "Hungarian.h"
#include "Metrics.h"

namespace HungarianWrapper{


    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations(
            const Eigen::PermutationMatrix<Eigen::Dynamic>& p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic>& p2) {

        long n1 = p1.size(), n2 = p2.size();
        Eigen::VectorXi combined(n1+n2);

        combined.segment(0,n1) = p1.indices().base();
        combined.segment(n1,n2) = (p2.indices().base().array()+n1);

        return Eigen::PermutationMatrix<Eigen::Dynamic>(combined);
    };

    Eigen::PermutationMatrix<Eigen::Dynamic> spinSpecificHungarian(
            const Reference &lhs,
            const Reference &rhs,
            bool flipSpinsQ = false) {
        assert(lhs.maximum_.typesVector() == rhs.maximum_.typesVector()
               && "The typesvectors must be identical as we assume ordered alpha and beta electrons.");
        assert(lhs.maximum_.positionsVector().numberOfEntities() == rhs.maximum_.positionsVector().numberOfEntities()
               && "The number of positions must be identical.");

        auto lhsCopy = lhs.maximum_.positionsVector();//TODO eliminate redundant copy by const_cast .slice
        //https://stackoverflow.com/questions/5008541/how-to-call-a-non-const-function-within-a-const-function-c
        auto rhsCopy = rhs.maximum_.positionsVector();

        auto nAlpha = lhs.maximum_.typesVector().countOccurence(Spin::alpha);
        auto nBeta = lhs.maximum_.typesVector().countOccurence(Spin::beta);

        Interval alphaElectrons{},betaElectrons{};
        if(flipSpinsQ){
            alphaElectrons = {nAlpha, nBeta};
            betaElectrons = {0, nAlpha};
        } else {
            alphaElectrons = {0, nAlpha};
            betaElectrons = {nAlpha, nBeta};
        }

        auto costMatrixAlpha = Metrics::positionalDistances(
                PositionsVector(lhsCopy.slice(alphaElectrons).dataRef()),
                PositionsVector(rhsCopy.slice(alphaElectrons).dataRef()));
        auto costMatrixBeta = Metrics::positionalDistances(
                PositionsVector(lhsCopy.slice(betaElectrons).dataRef()),
                PositionsVector(rhsCopy.slice(betaElectrons).dataRef()));

        auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);
        auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

        return combinePermutations(bestMatchAlpha,bestMatchBeta);
    };


};

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::set<Reference> references;
    std::vector<Sample> samples;

    RawDataReader reader(references,samples);
    reader.read("raw.bin");

    //for (auto ref : references) { ;
    //    console->info("{} {:03.6f}",ref.id_, ref.negLogSqrdProbabilityDensity_);
    //}

    double increment, distThresh = 0.01;
    if(!references.empty()) {
        increment = (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-5;
    } else {
        console->error("References are empty.");
        return false;
    }
    console->info("increment: {}",increment);
    console->flush();


    // first range


    auto uit = references.begin();
    auto lit = uit;
    while (lit != references.end()){

        lit = uit;
        uit = references.upper_bound(Reference((*lit).negLogSqrdProbabilityDensity_+increment));

        //TODO DOES THIS WORK FOR LIT == UIT? WHAT HAPPENS FOR UIT = LIT+1?

        console->info("{} ------ {}",(*lit).id_, (*lit).negLogSqrdProbabilityDensity_);
        console->info("{} ------ {}",(*uit).id_, (*uit).negLogSqrdProbabilityDensity_);

        auto it = lit;
        while(it++ != uit){
            // check hungarian
            auto bestMatch = HungarianWrapper::spinSpecificHungarian(*lit,*it);
            auto bestMatchFlip = HungarianWrapper::spinSpecificHungarian(*lit,*it,true);
            auto currentMaxCopy = (*it).maximum_;
            currentMaxCopy.permute(bestMatch);
            double mostDeviatingParticleDistance = Metrics::positionDistancesVector(
                    (*lit).maximum_.positionsVector(),
                    currentMaxCopy.positionsVector()).lpNorm<Eigen::Infinity>();

            currentMaxCopy = (*it).maximum_;
            currentMaxCopy.permute(bestMatch);
            double mostDeviatingParticleDistanceFlip = Metrics::positionDistancesVector(
                    (*lit).maximum_.positionsVector(),
                    currentMaxCopy.positionsVector()).lpNorm<Eigen::Infinity>();


            if((mostDeviatingParticleDistance <= distThresh)
            || (mostDeviatingParticleDistanceFlip <= distThresh) ){

                if(mostDeviatingParticleDistance <= mostDeviatingParticleDistanceFlip)
                    samples[(*it).id_].sample_.permute(bestMatch);
                else
                    samples[(*it).id_].sample_.permute(bestMatchFlip);

                // refs are identical
                (*lit).addAssociation((*it).id_);
                references.erase((*it));
                uit--;
            }
            uit++;
        }


    }
    console->flush();

}