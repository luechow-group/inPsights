//
// Created by Michael Heuer on 05.09.18.
//

#ifndef INPSIGHTS_LOGGER_H
#define INPSIGHTS_LOGGER_H

#include <string>
#include <spdlog/spdlog.h>

namespace Logger{

    const std::string name = "console";

    void initialize();
    std::shared_ptr<spdlog::logger> get();
}

#endif //INPSIGHTS_LOGGER_H
