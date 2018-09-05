//
// Created by Michael Heuer on 05.09.18.
//

#ifndef AMOLQCPP_LOGGER_H
#define AMOLQCPP_LOGGER_H

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace Logger{

    const std::string name = "console";

    void initialize();
}

#endif //AMOLQCPP_LOGGER_H
