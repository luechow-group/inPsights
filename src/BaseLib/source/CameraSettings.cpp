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

#include "CameraSettings.h"
#include <spdlog/spdlog.h>
#include <string>

namespace Settings {
    Camera::Camera()
            : ISettings(VARNAME(Camera)) {
        zoom.onChange_.connect(
                [&](int value) {
                    if (not (value >= 0 and value < 999))
                        throw std::invalid_argument("The zoom must be in bewteen ["
                        + std::to_string(0) + "," + std::to_string(999) + "]");
                });
        pan.onChange_.connect(
                [&](int value) {
                    int limit = 180;
                    if (std::abs(value) > limit)
                        throw std::invalid_argument("The pan angle must be inbetween ["
                        + std::to_string(-limit) + "," + std::to_string(limit) + "]");
                });
        tilt.onChange_.connect(
                [&](int value) {
                    int limit = 180;
                    if (std::abs(value) > limit)
                        throw std::invalid_argument("The tilt angle must be inbetween ["
                                                    + std::to_string(-limit) + "," + std::to_string(limit) + "]");
                });
        roll.onChange_.connect(
                [&](int value) {
                    int limit = 180;
                    if (std::abs(value) > limit)
                        throw std::invalid_argument("The roll angle must be inbetween ["
                                                    + std::to_string(-limit) + "," + std::to_string(limit) + "]");
                });
    }

    Camera::Camera(const YAML::Node &node)
            : Camera() {
        intProperty::decode(node[className], zoom);
        intProperty::decode(node[className], pan);
        intProperty ::decode(node[className], tilt);
        intProperty ::decode(node[className], roll);
    }

    void Camera::appendToNode(YAML::Node &node) const {
        node[className][zoom.name()] = zoom.get();
        node[className][pan.name()] = pan.get();
        node[className][tilt.name()] = tilt.get();
        node[className][roll.name()] = roll.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::Camera)

namespace Camera {
    Settings::Camera settings = Settings::Camera();
}