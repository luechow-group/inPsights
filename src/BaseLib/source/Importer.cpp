// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later


#include <numeric>
#include "Importer.h"
#include <iostream>
#include <exception>

Importer::Importer(const std::string &filename)
        : lines_(read_lines(filename)) {}

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

std::string Importer::strip(const std::string &s) const {
    auto vec = split(s);
    return std::accumulate(vec.begin(), vec.end(), std::string(""));
}

std::string Importer::strip(const std::string &s, char delim) const {

    auto vec = split(s,delim);
    for (auto it = vec.begin(); it != vec.end(); it++){
        if( *it == std::to_string(delim)) vec.erase(it);
    }
    return std::accumulate(vec.begin(), vec.end(), std::string(""));
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

//template<class OutIt>
//void Importer::read_lines(std::istream& is, OutIt dest)
//{
//    typedef std::istream_iterator<Line> InIt;
//    std::copy(InIt(is), InIt(), dest);
//}
std::vector<std::string> Importer::read_lines(const std::string& filename) {

    std::ifstream file;
    try {
        file.open(filename);
        if(!file.is_open()){
            std::cerr << "File could not be be found" << std::endl;
            throw std::exception();
        }

        while (!file.eof()) file.get();
        file.close();
    }
    catch (std::ifstream::failure& e) {
        std::cerr << "Exception opening/reading/closing file\n";
    }

    file.open (filename);
    std::vector<std::string> lines;
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            lines.push_back(line);
        }
    };
    file.close();

    return lines;
}



