/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
                             "  " + inPsights::version + "\n"
                             "\n"
                             " Main author:\n"
                             "  Michael A. Heuer, RWTH Aachen University\n"
                             "\n"
                             " With contributions from:\n"
                             "  Leonard Reuter\n";
}
#endif //INPSIGHTS_LOGO_H
