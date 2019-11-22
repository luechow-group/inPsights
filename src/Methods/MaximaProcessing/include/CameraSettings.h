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

#ifndef INPSIGHTS_CAMERASETTINGS_H
#define INPSIGHTS_CAMERASETTINGS_H

#include <ISettings.h>

namespace Settings {
    class Camera : public ISettings {
    public:
        Property<float> distance = {8.0f, VARNAME(distance)};
        Property<float> pan = {0.0f, VARNAME(pan)};
        Property<float> tilt = {45.0f, VARNAME(tilt)};
        Property<float> roll = {0.0f, VARNAME(roll)};

        Camera();
        explicit Camera(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::Camera)

namespace Camera {
    extern Settings::Camera settings;
}

#endif //INPSIGHTS_CAMERASETTINGS_H
