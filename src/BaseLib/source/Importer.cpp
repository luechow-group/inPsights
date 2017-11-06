//
// Created by Michael Heuer on 02.11.17.
//

#include <iterator>
#include <vector>
#include "Importer.h"
#include "ElementInfo.h"



Importer::Importer(const std::string &filename)
        : filename_(filename) {
    file_.open(filename);
    read_lines(file_, std::back_inserter(lines_));
}

Importer::~Importer() {
    file_.close();
}

std::string Importer::getLine(unsigned long idx) const {
    return lines_[idx];
}

// splits whitespace and ignores repeated whitespaces
std::vector<std::string> Importer::split(const std::string &s) const {
    std::istringstream iss(s);
    return std::vector<std::string>({std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>{}});
}


// splits any delimiter type
template<typename Out>
void Importer::split(const std::string &s, char delimiter, Out result) const {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delimiter)) {
        *(result++) = item;
    }
}

std::vector<std::string> Importer::split(const std::string &s, char delim) const {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

class Line {
    std::string data;
public:
    friend std::istream &operator>>(std::istream &is, Line &l) {
        std::getline(is, l.data);
        return is;
    }
    operator std::string() const { return data; }
};

template<class OutIt>
void Importer::read_lines(std::istream& is, OutIt dest)
{
    typedef std::istream_iterator<Line> InIt;
    std::copy(InIt(is), InIt(), dest);
}

AmolqcImporter::AmolqcImporter(const std::string &filename)
        : Importer(filename) {}


RefFileImporter::RefFileImporter(const std::string &filename)
        : AmolqcImporter(filename) {

    std::cout << getLine(0) << std::endl;

    numberOfNuclei_ = std::stoul(split(getLine(0))[0]);
    numberOfElectrons_ = std::stoul(split(getLine(numberOfNuclei_+2))[10]);
    numberOfAlphaElectrons_ = std::stoul(split(getLine(numberOfNuclei_+2))[8]);

            std::vector<std::string> summaryLineElements = split(getLine(numberOfNuclei_+1));

    numberOfSuperstructures_ = std::stoul(summaryLineElements[0]);
    totalNumberOfMaxima_ = std::stoul(summaryLineElements[2]);
    maximalNumberOfSubstructures = std::stoul(summaryLineElements[3]);
    
    countSubstructures();
}

AtomCollection RefFileImporter::getAtomCollection() {
    AtomCollection atomCollection;

    for (unsigned i = 1; i <= numberOfNuclei_; ++i) {
        std::vector<std::string> lineElements = split(getLine(i));
        Elements::ElementType elementType = Elements::ElementInfo::elementTypeForSymbol(lineElements[1]);
        double x = std::stod(lineElements[2]);
        double y = std::stod(lineElements[3]);
        double z = std::stod(lineElements[4]);
        atomCollection.addAtom(x,y,z,elementType);
    }

    return atomCollection;
}

void RefFileImporter::countSubstructures() {
    unsigned long currentLineIdx = numberOfNuclei_+2;
    unsigned long next = numberOfElectrons_+2;

    unsigned long sumOfMaximaNumbersTillCurrent = 0;
    unsigned long sumOfMaximaNumbersWithCurrent = 0;

    std::string currentLine = getLine(currentLineIdx);
    unsigned long firstLineOfSuperstructure = currentLineIdx;

    unsigned long k = 0;
    unsigned long m = 0;
    unsigned long k_last = 1;
    unsigned long m_last = 1;

    while (!currentLine.empty()){

        std::vector<std::string> currentLineElements = split(currentLine);
        k = std::stoul(currentLineElements[1]);
        m = std::stoul(currentLineElements[2]);

        // if new superstructure reached
        if ( k > k_last){
            sumOfMaximaNumbersWithCurrent = std::stoul(currentLineElements[6]);

            substructuresData_.emplace_back(
                    SubstructureDataEntry(firstLineOfSuperstructure, m_last, sumOfMaximaNumbersTillCurrent));
            firstLineOfSuperstructure = currentLineIdx;
            k_last = k;
            m_last = 1;
        }
            // else currentLine contains another substructure
        else {
            sumOfMaximaNumbersWithCurrent = sumOfMaximaNumbersTillCurrent;
            sumOfMaximaNumbersWithCurrent += std::stoul(currentLineElements[6]);
            m_last = m;
        };

        currentLineIdx +=  next;
        currentLine = getLine(currentLineIdx);
        sumOfMaximaNumbersTillCurrent = sumOfMaximaNumbersWithCurrent;
     }

    // add last superstructure
    if (k > 0){
        m_last;
        substructuresData_.emplace_back(
                SubstructureDataEntry(firstLineOfSuperstructure, m_last, sumOfMaximaNumbersTillCurrent));
    }
}

unsigned long RefFileImporter::calculateLine(unsigned long k, unsigned long m) const {
    assert( k > 0 && "k value must be greater than zero");
    assert( m > 0 && "m value must be greater than zero");
    assert( k <= numberOfSuperstructures_ && "k value must be smaller than kmax");
    assert( m <= substructuresData_[k].startingLine_  && "m value must be smaller than mmax");

    unsigned long start = substructuresData_[k-1].startingLine_;
    unsigned long linesToSkip = (m-1)*(numberOfElectrons_+2);
    return start + linesToSkip;
}

SpinTypeCollection RefFileImporter::getSpinTypeCollection() const {
    SpinTypeCollection spinTypeCollection;
    for (unsigned long i = 0; i < numberOfElectrons_; ++i) {
        Spin::SpinType spinType;
        if (i < numberOfAlphaElectrons_) spinType = Spin::SpinType::alpha;
        else  spinType = Spin::SpinType::beta;

        spinTypeCollection.append(spinType);
    }
    return spinTypeCollection;
}

ParticleCollection RefFileImporter::getParticleCollection(unsigned long k, unsigned long m) const {

    unsigned long startLine = calculateLine(k,m)+2;
    ParticleCollection particleCollection;

    for (unsigned long i = 0; i < numberOfElectrons_; ++i) {
        std::vector<std::string> lineElements = split(getLine(startLine+i));
        double x = std::stod(lineElements[0]);
        double y = std::stod(lineElements[1]);
        double z = std::stod(lineElements[2]);

        particleCollection.append(Particle(x,y,z));
    }
    return particleCollection;
}

//maxref
unsigned long RefFileImporter::getNumberOfMaxima(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+2;
    std::vector<std::string> lineElements = split(getLine(startLine));
    return std::stoul(lineElements[6]);
}

//maxref
double RefFileImporter::getNegativeLogarithmizedProbabilityDensity(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+2;
    std::vector<std::string> lineElements = split(getLine(startLine));
    return std::stod(lineElements[4]);
}

ElectronCollection RefFileImporter::getElectronCollection(unsigned long k, unsigned long m) const {
    return ElectronCollection(this->getParticleCollection(k,m), this->getSpinTypeCollection());
}

ElectronCollections RefFileImporter::getElectronCollections(unsigned long k) const {
    unsigned long numberOfSubstructures = substructuresData_[k].numberOfSubstructures_;

    std::vector<ParticleCollection> particleCollectionVector;
    for (int m = 1; m <= numberOfSubstructures; ++m) {
        particleCollectionVector.emplace_back(this->getParticleCollection(k,m));
    }
    return ElectronCollections(particleCollectionVector,this->getSpinTypeCollection());

}

