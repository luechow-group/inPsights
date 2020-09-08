// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ERRORHANDLING_H
#define INPSIGHTS_ERRORHANDLING_H

#include <stdexcept>

class NotImplemented : public std::logic_error
{
public:
    NotImplemented();
};

#endif //INPSIGHTS_ERRORHANDLING_H
