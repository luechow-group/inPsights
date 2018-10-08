//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_RAWDATAREADER_H
#define AMOLQCPP_RAWDATAREADER_H

#include <BinaryFileReader.h>
#include "Reference.h"
#include "Sample.h"

class RawDataReader : public BinaryFileReader{
public:
    explicit RawDataReader(
            std::vector<Reference>& references,
            std::vector<Sample>& samples,
            int recordDelimiterLength = 4);

    void read(const std::string& fileName) override;

private:
    std::vector<Reference>& references_;
    std::vector<Sample>& samples_;
};

#endif //AMOLQCPP_RAWDATAREADER_H