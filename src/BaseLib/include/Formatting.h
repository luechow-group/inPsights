// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_FORMATTING_H
#define INPSIGHTS_FORMATTING_H

#include <spdlog/fmt/bundled/format.h>
#include <Eigen/Core>
#include <cstdint>

// this user-defined type formatter is needed for settings of VoxelCubeGeneration
// it is written based on a documentation example: https://fmt.dev/latest/api.html#udt
// it is needed, because settings defaults are printed by speedlog in Property.h

template<typename T>
struct fmt::formatter<Eigen::Matrix<T,3,1>> {
    // Dummy parse routine
    constexpr auto parse(format_parse_context& ctx) {
        return ctx.end();
    }

    template <typename FormatContext>
    auto format(const Eigen::Matrix<T,3,1>& matrix, FormatContext& ctx) {
        return format_to(
                ctx.out(),
                "[{}, {}, {}]", matrix[0], matrix[1], matrix[2]);
    }
};

#endif //INPSIGHTS_FORMATTING_H
