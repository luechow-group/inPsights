//
// Created by Michael Heuer on 05.09.18.
//

#include <Logger.h>

namespace Logger{
    void initialize(){
        auto console = spdlog::stdout_color_st(name);
        console->info("Welcome to spdlog!");
    };

}