// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BuildInfo.h>

std::string inPsights::version() {
    std::string version;
#ifdef INPSIGHTS_VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
    version = TOSTRING(INPSIGHTS_VERSION_INFO);
    // remove quotation marks
    version.substr(1, version.size() - 2);
#else
    version = "unkown!";
#endif
    return version;
}

std::string inPsights::compiler() {
    std::string compiler;
#ifdef INPSIGHTS_CXX_COMPILER_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
    compiler = TOSTRING(INPSIGHTS_CXX_COMPILER_INFO);
    // remove quotation marks
    compiler.substr(1, compiler.size() - 2);
#else
    compiler = "unkown!";
#endif
    return compiler;
}