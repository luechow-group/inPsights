// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ENUMERATE_H
#define INPSIGHTS_ENUMERATE_H

#include <tuple>

template <typename T,
        typename TIter = decltype(std::begin(std::declval<T>())),
        typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable)
{
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator != (const iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper
    {
        T iterable;
        auto begin() { return iterator{ 0, std::begin(iterable) }; }
        auto end() { return iterator{ 0, std::end(iterable) }; }
    };
    return iterable_wrapper{ std::forward<T>(iterable) };
}

#endif //INPSIGHTS_ENUMERATE_H
