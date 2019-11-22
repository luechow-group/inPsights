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

#include "CameraSettings.h"
#include <spdlog/spdlog.h>
#include <string>

namespace Settings {
    Camera::Camera()
            : ISettings(VARNAME(Camera)) {
        distance.onChange_.connect(
                [&](float value) {
                    if (not (value >= 0.0f and value < std::numeric_limits<float>::max()))
                        throw std::invalid_argument("The camera distance can not be negative.");
                });
        pan.onChange_.connect(
                [&](float value) {
                    float limit = 180.0f;
                    if (std::abs(value) > limit)
                        throw std::invalid_argument("The pan angle must be inbetween ["
                        + std::to_string(-limit) + "," + std::to_string(limit) + "]");
                });
        tilt.onChange_.connect(
                [&](float value) {
                    float limit = 90.0f;
                    if (std::abs(value) > limit)
                        throw std::invalid_argument("The tilt angle must be inbetween ["
                                                    + std::to_string(-limit) + "," + std::to_string(limit) + "]");
                });
        roll.onChange_.connect(
                [&](int value) {
                    float limit = 180.0f;
                    if (std::abs(value) > limit)
                        throw std::invalid_argument("The roll angle must be inbetween ["
                                                    + std::to_string(-limit) + "," + std::to_string(limit) + "]");
                });
    }

    Camera::Camera(const YAML::Node &node)
            : Camera() {
        floatProperty::decode(node[className], distance);
        floatProperty::decode(node[className], pan);
        floatProperty ::decode(node[className], tilt);
        floatProperty ::decode(node[className], roll);
    }

    void Camera::appendToNode(YAML::Node &node) const {
        node[className][distance.name()] = distance.get();
        node[className][pan.name()] = pan.get();
        node[className][tilt.name()] = tilt.get();
        node[className][roll.name()] = roll.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::Camera)

namespace Camera {
    Settings::Camera settings = Settings::Camera();
}