//
// Created by Michael Heuer on 2018-12-31.
//

#include <Logger.h>

namespace Logger{
    std::shared_ptr<spdlog::logger> console = spdlog::stdout_color_st("console");
}