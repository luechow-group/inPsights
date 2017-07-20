//
// Created by Morian Sonneton 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_INPUTOUTPUT_H
#define LOCALSPINMULTIPLICITY_INPUTOUTPUT_H
#include <string>
#include <vector>
#include <fstream>

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
