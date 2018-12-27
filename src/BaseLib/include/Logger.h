//
// Created by Michael Heuer on 05.09.18.
//

#ifndef INPSIGHTS_LOGGER_H
#define INPSIGHTS_LOGGER_H

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace Logger{
    inline std::shared_ptr<spdlog::logger> console = {spdlog::stdout_color_st("console")};
}

#endif //INPSIGHTS_LOGGER_H
