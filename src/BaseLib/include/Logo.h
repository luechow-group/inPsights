// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_LOGO_H
#define INPSIGHTS_LOGO_H

#include <Version.h>

namespace inPsights {
    const std::string logo = "                  ▀▓▓▓▓▀ ▄▄■■■■▄▄\n"
                             "                    ▓▓  █ ▄▄▄▄▄  █\n"
                             "▐▌         ▀▓▓▌     ▓▓ █  ▓▓▓▀    █        ▐▌\n"
                             "            ▐▓▓     ▓▓  █ ▓▓▓    █         ▐▌         ▐▌\n"
                             "▐▌  ▐▌▄▀▀▀▄  ▓▓     ▓▓   ▀▀■■■■▀▀ ▄▀▀▀▄▐▌  ▐▌▄▀▀▀▄  ▀▀▐▌▀▀  ▄▀▀▀▄\n"
                             "▐▌  ▐▌    ▐▌ ▀▓▓    ▓▓    ▓▓▀ ░░ ▐▌    ▐▌  ▐▌    ▐▌   ▐▌   ▐▌\n"
                             "▐▌  ▐▌    ▐▌   ▀▓   ▓▓   ▓▀    ░░▐▌    ▐▌  ▐▌    ▐▌   ▐▌    ▀■■■▄\n"
                             "▐▌  ▐▌    ▐▌     ▀▓▄▓▓▄▓▀        ▐▌    ▐▌  ▐▌    ▐▌   ▐▌        ▐▌\n"
                             "▐▌  ▐▌    ▐▌        ▓▓            ▀▄▄▄▀▐▌  ▐▌    ▐▌    ■■   ▀▄▄▄▀\n"
                             "                  ▄▓▓▓▓▄               ▐▌\n"
                             "                                  ▀■■■■▀\n"
                             "Version:\n"
                             "  " + inPsights::version() + "\n"
                             "Compiler:\n"
                             "  " + inPsights::compiler() + "\n"
                             "\n"
                             " Main author:\n"
                             "  Michael A. Heuer, RWTH Aachen University\n"
                             "\n"
                             " With contributions from:\n"
                             "  Leonard Reuter, RWTH Aachen University\n";
}
#endif //INPSIGHTS_LOGO_H
