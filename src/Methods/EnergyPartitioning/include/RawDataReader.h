//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_RAWDATAREADER_H
#define AMOLQCPP_RAWDATAREADER_H

#include <BinaryFileReader.h>

class RawDataReader : public BinaryFileReader{
public:
    explicit RawDataReader(int recordDelimiterLength = 4);

    bool read(const std::string& fileName) override;
};

#endif //AMOLQCPP_RAWDATAREADER_H
