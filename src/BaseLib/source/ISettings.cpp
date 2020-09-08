// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <ISettings.h>

Settings::ISettings::ISettings(const std::string& className)
: className(className){}

std::string Settings::ISettings::name() const {
    return className;
}
