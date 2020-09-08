// Copyright (C) 2019 Michael Heuer.
// Copyright (C) 2020 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

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
