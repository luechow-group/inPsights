// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_IMPORTER_H
#define INPSIGHTS_IMPORTER_H

#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <sstream>

class StringManipulator{
public:
    std::vector<std::string> split(const std::string &s){
        std::istringstream iss(s);
        return std::vector<std::string>({std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>{}});
    }
};

class Importer{

public:
    explicit Importer(const std::string& filename);

    std::string getLine(unsigned long idx) const;

    // TODO move to a string helper class?
    template<typename Out>
    void split(const std::string &s, char delimiter, Out result) const;

    std::vector<std::string> split(const std::string &s) const;
    std::vector<std::string> split(const std::string &s, char delim) const;

    std::string strip(const std::string &s) const;
    std::string strip(const std::string &s, char delim) const;

private:
    std::vector<std::string> read_lines(const std::string& file);

protected:
    std::vector<std::string> lines_;
};

class SubstructureDataEntry{
public:
    SubstructureDataEntry(unsigned long startingLine,
                          unsigned long numberOfSubstructures,
                          unsigned long totalNumberOfMaxima = 1)
            : startingLine_(startingLine),
              numberOfSubstructures_(numberOfSubstructures),
              totalNumberOfMaxima_(totalNumberOfMaxima) {};

    unsigned long startingLine_, numberOfSubstructures_, totalNumberOfMaxima_;
};

#endif //INPSIGHTS_IMPORTER_H
