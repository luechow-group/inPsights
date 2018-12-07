//
// Created by Michael Heuer on 05.09.18.
//

#include <Logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace Logger{
    void initialize(){
        auto console = spdlog::stdout_color_st(name);
    };

    std::shared_ptr<spdlog::logger> get() {
        return spdlog::get(Logger::name);
    }
}