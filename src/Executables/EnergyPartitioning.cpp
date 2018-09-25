//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Sample.h>
#include <Logger.h>
#include <HungarianHelper.h>
#include <algorithm>
#include <DensityBasedScan.h>
#include <CoulombPotential.h>
#include <Statistics.h>

#include <MoleculeWidget.h>
#include <ElectronicWaveFunction.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <QApplication>


class GlobalIdentiySorter{
public:

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double increment, double distThresh)
    :
    references_(references),
    samples_(samples),
    increment_(increment),
    distThresh_(distThresh),
    console(spdlog::get(Logger::name))
    {}

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double distThresh = 0.01)
    :
    GlobalIdentiySorter(references, samples, (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-7, distThresh)
    {}


    bool sort(){
        if(references_.empty()) {
            console->error("References are empty.");
            return false;
        }
        else if (references_.size() == 1) {
            console->warn("No sorting because only one reference was found.");
            return true; // no sorting
        }
        std::sort(references_.begin(),references_.end());

        auto lit = references_.begin();
        auto uit = references_.begin();

        while (lit != references_.end()){
            uit = std::upper_bound(lit,references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            //uit = references_.upper_bound(Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            auto it = lit;
            it++; // start with the element next to lit

            while(it != uit) subLoop(lit,it);

            lit = uit;
        }

        return true;
    }

private:
    void subLoop(std::vector<Reference>::iterator& lit, std::vector<Reference>::iterator& it){
        //TODO maybe calculate only alpha electron distance and skip beta electron hungarian if dist is too high already

        auto bestMatch = Metrics::spinSpecificBestMatchNorm((*it).maximum_, (*lit).maximum_);


        if((*lit).maximum_.typesVector().multiplicity() == 1) { // consider spin flip
            auto bestMatchFlipped = Metrics::spinSpecificBestMatchNorm((*it).maximum_, (*lit).maximum_, true);

            console->info("{},{}",bestMatch.first,bestMatchFlipped.first);
            if( (bestMatch.first <= distThresh_) || (bestMatchFlipped.first <= distThresh_) ){
                console->info("MERGE");

                if(bestMatch.first <= bestMatchFlipped.first)
                    addReference(lit, it, bestMatch.second);
                else
                    addReference(lit, it, bestMatchFlipped.second);

            } else
                it++;
        }
        else {  // don't consider spin flip
            if( (bestMatch.first <= distThresh_) )
                addReference(lit, it, bestMatch.second);
            else
                it++;
        }
    }

    //TODO find better name
    void addReference(const std::vector<Reference, std::allocator<Reference>>::iterator &lit,
                      std::vector<Reference, std::allocator<Reference>>::iterator &it,
                      const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {
        samples_[(*it).id_].sample_.permute(bestMatch);

        (*lit).addAssociation((*it).id_);
        (*lit).associatedSampleIds_.insert(
                        (*lit).associatedSampleIds_.end(),
                        make_move_iterator((*it).associatedSampleIds_.begin()),
                        make_move_iterator((*it).associatedSampleIds_.end())
                );

        it = references_.erase(it); // returns the iterator of the following element
    }

    std::vector<Reference>& references_;
    std::vector<Sample>& samples_;
    double increment_,distThresh_;
    std::shared_ptr<spdlog::logger> console;

};


class GlobalSimilaritySorter {
public:

    GlobalSimilaritySorter(
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double increment, double distThresh)
    :
    references_(references),
    similarReferencesVector_(similarReferencesVector),
    increment_(increment),
    distThresh_(distThresh),
    console(spdlog::get(Logger::name))
    {}

    GlobalSimilaritySorter(
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double distThresh = 0.1)
    :
    GlobalSimilaritySorter(references, similarReferencesVector,
            (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-5, distThresh)
    {}


    bool sort(){

        if(similarReferencesVector_.empty())
            similarReferencesVector_.emplace_back(SimilarReferences(references_.begin()));

        // for all references in range
        // for (auto refIt = references_.begin()+1; refIt != references_.end(); ++refIt) {
        auto lit = references_.begin();
        auto uit = references_.begin();
        uit = std::upper_bound(references_.begin(),references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));


        while (lit != references_.end()) {
            uit = std::upper_bound(lit,references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));

            for (auto it = lit; it != uit;  ++it ) {

                bool isSimilarQ = false;
                for (auto &simRefs : similarReferencesVector_) {
                    // check if refIt is similar to simRefs representative reference

                    // calc perm r->sr so that we can store them in similar reference
                    auto bestMatch = Metrics::bestMatchNorm<Eigen::Infinity>(
                            (*it).maximum_.positionsVector(),
                            (*simRefs.representativeReferenceIterator).maximum_.positionsVector());

                    //console->info("{}", dist);
                    if (bestMatch.first < distThresh_) {
                        simRefs.similarReferences_.emplace_back(SimilarReference(it, bestMatch.second));

                        isSimilarQ = true;
                    }
                }
                if (!isSimilarQ) similarReferencesVector_.emplace_back(SimilarReferences(it));
            }
            lit = uit;
        }
        return true;
    }
private:

    std::vector<Reference>& references_;
    std::vector<SimilarReferences>& similarReferencesVector_;
    double increment_,distThresh_;
    std::shared_ptr<spdlog::logger> console;
};


int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);
    double identicalDistThresh = 0.01;
    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;

    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("Ethane.bin");
    console->info("number of refs {}",globallyIdenticalMaxima.size());
    auto numberOfElectrons = globallyIdenticalMaxima[0].maximum_.numberOfEntities();


    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.sort();
    console->info("finished id sort");
    console->flush();


    console->info("total elems {}",std::distance(globallyIdenticalMaxima.begin(),globallyIdenticalMaxima.end()));


    double similarDistThresh = 0.1;
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(globallyIdenticalMaxima, similarReferencesVector,similarDistThresh);

    globalSimilaritySorter.sort();

    console->info("total elems {}",similarReferencesVector.size());


    std::vector<SimilarReference> clusteredReferences;

    /*Clustering*/
    // make data
    std::vector<Eigen::VectorXd> data;
    for(auto& simRefVector : similarReferencesVector){
        data.push_back((*simRefVector.representativeReferenceIterator).maximum_.positionsVector().asEigenVector());
    }

    // needs bestMatchDistance
    DensityBasedScan<double> dbscan(data);

    auto nClusters = dbscan.findClusters(0.20001, 5);
    auto result = dbscan.getLabels();


    /*Energy Partitioning*/
    Statistics::RunningStatistics<Eigen::VectorXd> EkinStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> EpotStats;


    Eigen::VectorXd ekin;
    Eigen::MatrixXd epot;
    size_t sampleId, count, simCount, totalCount;

    totalCount = 0;

    int i=0;
    for(auto& simRefVector : similarReferencesVector){
        console->info("{}: contained similar refs {}",i,simRefVector.similarReferences_.size());
        i++;
        simCount = 0;

        EkinStats.reset();
        EpotStats.reset();

        // Representative reference
        count = 1+(*simRefVector.representativeReferenceIterator).associatedSampleIds_.size();
        simCount += count;
        sampleId = (*simRefVector.representativeReferenceIterator).id_;


        ekin = samples[sampleId].kineticEnergies_;
        epot = CoulombPotential::energies(samples[sampleId].sample_);

        console->info("rep ref count {}" ,count);
        EkinStats.add(ekin,unsigned(count));
        EpotStats.add(epot,unsigned(count));

        // Iterate over references being similar to the representative reference.
        for(const auto& simRef : simRefVector.similarReferences_){


            //TODO CONTINUE HERE
            globallyIdenticalMaxima[0];

            // Hier werden irgendwo maxima counts vergessen#
            count = 1+(*simRef.it_).associatedSampleIds_.size();
            simCount += count;
            sampleId = (*simRef.it_).id_;

            auto sampleCopy = samples[sampleId].sample_;
            sampleCopy.permute(simRef.perm_);

            ekin = simRef.perm_ * (samples[sampleId].kineticEnergies_);
            epot = CoulombPotential::energies(sampleCopy);

            EkinStats.add(ekin,unsigned(count));
            EpotStats.add(epot,unsigned(count));
        }
        console->info("sim ref count {}",simCount);
        std::cout <<"mean: ("<<EkinStats.getWeightedSum()<<")\n"<< EkinStats.mean().transpose() << std::endl;
        if(simCount >=2)
            std::cout <<"stdv:\n"<< EkinStats.standardDeviation().transpose() << std::endl<< std::endl;

        totalCount += simCount;
    }
    console->info("overall count {}",totalCount);




    std::string wavefunctionFilename = "Ethane-em.wf";
    auto wf = ElectronicWaveFunction::getInstance(wavefunctionFilename);
    auto atoms = wf.getAtomsVector();

    // Visualization
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomsVector3D(root, atoms);

    auto ev1 = (*similarReferencesVector.at(1).representativeReferenceIterator).maximum_;
    auto perm = similarReferencesVector.at(1).similarReferences_.at(0).perm_;
    auto ev2 = (*similarReferencesVector.at(1).similarReferences_.at(0).it_).maximum_;
    ev2.permute(perm);

    ElectronsVector3D(root, atoms, ev1, true);
    ElectronsVector3D(root, atoms, ev2, true);

    return app.exec();


//TODO CHECK PERMUTATION
// CHECK ADDITIONAL +1
//TODO DBSCAN

};
