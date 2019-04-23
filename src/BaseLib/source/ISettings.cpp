//
// Created by Michael Heuer on 2019-04-23.
//

#include <ISettings.h>

Settings::ISettings::ISettings(const std::string& className)
: className(className){}

std::string Settings::ISettings::getClassname() const {
    return className;
}
