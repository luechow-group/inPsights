//
// Created by Michael Heuer on 05.09.18.
//

#ifndef INPSIGHTS_LOGGER_H
#define INPSIGHTS_LOGGER_H

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace Logger{
    extern std::shared_ptr<spdlog::logger> console;
}

#endif //INPSIGHTS_LOGGER_H
