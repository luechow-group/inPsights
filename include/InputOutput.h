//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_INPUTOUTPUT_H
#define LOCALSPINMULTIPLICITY_INPUTOUTPUT_H
#include <string>
#include <vector>
#include <fstream>

/*
 * This class handles the used filestreams.
 * It is intended to also use this class for Output, which is not implemented yet.
 */
class InputOutput {
public:
    InputOutput(const std::string &filename, bool isRead);
    virtual ~InputOutput();
    void openFile(const std::string &filename, bool isRead);
    void closeAllFiles();
protected:
    std::vector<std::fstream> streams;
private:
    std::vector<std::string> filenames;
};


#endif //LOCALSPINMULTIPLICITY_INPUTOUTPUT_H
