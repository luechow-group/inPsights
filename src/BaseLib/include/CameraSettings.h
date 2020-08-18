/* Copyright (C) 2019 Michael Heuer.
 * Copyright (C) 2020 Leonard Reuter.
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
        Property<int> zoom = {100, VARNAME(zoom)};
        Property<int> pan = {0, VARNAME(pan)};
        Property<int> tilt = {0, VARNAME(tilt)};
        Property<int> roll = {0, VARNAME(roll)};

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
