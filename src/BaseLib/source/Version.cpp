// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Version.h>

std::string inPsights::version() {
    std::string version;
#ifdef INPSIGHTS_VERSION
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
    version = TOSTRING(INPSIGHTS_VERSION);

#else
    version = "unkown!";
#endif
    return version;
}
