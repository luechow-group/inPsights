//
// Created by Michael Heuer on 2018-12-27.
//

#ifndef INPSIGHTS_PROPERTY_H
#define INPSIGHTS_PROPERTY_H

#include <utility>
#include <boost/signals2/signal.hpp>
#include <iostream>
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>
#include <spdlog/spdlog.h>
// A Property is a encapsulates a value and may inform
// you on any changes applied to this value.

template<typename T>
class Property {

public:
    // Properties for built-in types are automatically
    // initialized to 0. See template spezialisations
    // at the bottom of this file
    Property(std::string name = "")
            : value_(0), connection_(nullptr), connectionId_(-1), name_(std::move(name)) {}

    Property(T const &val, std::string name = "")
            : value_(val), connection_(nullptr), connectionId_(-1), name_(std::move(name)) {}

    Property(T &&val, std::string name = "")
            : value_(std::move(val)), connection_(nullptr), connectionId_(-1), name_(std::move(name)) {}

    Property(Property<T> const &toCopy, std::string name = "")
            : value_(toCopy.value_), connection_(nullptr), connectionId_(-1), name_(std::move(name)) {}

    Property(Property<T> &&toCopy, std::string name = "")
            : value_(std::move(toCopy.value_)), connection_(nullptr), connectionId_(-1), name_(std::move(name)) {}


    // sets the Property to a new value. beforeChange() and
    // onChange() will be emitted.
    virtual void set(T const &value) {
        if (value != value_) {
            beforeChange_(value_);
            value_ = value;
            onChange_(value_);
        }
    }

    // sets the Property to a new value. beforeChange() and
    // onChange() will not be emitted
    void setWithNoEmit(T const &value) {
        value_ = value;
    }

    // emits beforeChange() and onChange() even if the value
    // did not change
    void touch() {
        beforeChange_(value_);
        onChange_(value_);
    }

    // returns the internal value
    virtual T const &get() const { return value_; }
    
    // if there are any Properties connected to this Property,
    // they won't be notified of any further changes
    virtual void disconnect_auditors() {
        onChange_.disconnect_all_slots();
        beforeChange_.disconnect_all_slots();
    }

    // assigns the value of another Property
    virtual Property<T> &operator=(Property<T> const &rhs) {
        set(rhs.value_);
        return *this;
    }

    // assigns a new value to this Property
    virtual Property<T> &operator=(T const &rhs) {
        set(rhs);
        return *this;
    }

    // compares the values of two Properties
    bool operator==(Property<T> const &rhs) const {
        return Property<T>::get() == rhs.get();
    }

    bool operator!=(Property<T> const &rhs) const {
        return Property<T>::get() != rhs.get();
    }

    // compares the values of the Property to another value
    bool operator==(T const &rhs) const { return Property<T>::get() == rhs; }

    bool operator!=(T const &rhs) const { return Property<T>::get() != rhs; }

    // returns the value of this Property
    T const &operator()() const { return Property<T>::get(); }

    std::string name() const { return name_; }

private:
    T value_;
    Property<T> const *connection_;
    int connectionId_;
    std::string name_;
public:
    boost::signals2::signal<void(T)> onChange_;
    boost::signals2::signal<void(T)>beforeChange_;
};

// stream operators
template<typename T>
std::ostream &operator<<(std::ostream &out_stream, Property<T> const &val) {
    out_stream << val.name() << " " << val.get();
    return out_stream;
}

template<typename T>
std::istream &operator>>(std::istream &in_stream, Property<T> &val) {
    std::string name;
    in_stream >> name;

    val.name(name);
    T tmp;
    in_stream >> tmp;
    val.set(tmp);
    return in_stream;
}


namespace YAML {
    template<typename T>
    struct convert<Property<T>> {
        static Node encode(const Property<T> &t) {
            Node node;
            node[t.name()] = t.get();
            return node;
        }

        static bool decode(const Node &node, Property<T> &rhs) {
            if(!node.IsMap() || !node[rhs.name()]) {
                spdlog::warn("Property \"{0}\" was not found. Using preset value: {1}", rhs.name(), rhs.get());
                return false;
            }

            rhs = node[rhs.name()].template as<T>();
            return true;
        }
    };

    template<typename T>
    Emitter &operator<<(Emitter &out, const Property<T> &t) {
        out << Key << t.name() << Value << t.get();
        return out;
    };
}

using doubleProperty        = YAML::convert<Property<double>>;
using floatProperty         = YAML::convert<Property<float>>;
using shortProperty         = YAML::convert<Property<short>>;
using intProperty           = YAML::convert<Property<int>>;
using longProperty          = YAML::convert<Property<long>>;
using charProperty          = YAML::convert<Property<char>>;
using stringProperty        = YAML::convert<Property<std::string>>;
using size_tProperty        = YAML::convert<Property<size_t>>;
using unsignedProperty      = YAML::convert<Property<unsigned>>;
using unsignedShortProperty = YAML::convert<Property<unsigned short>>;
using boolProperty          = YAML::convert<Property<bool>>;

#endif //INPSIGHTS_PROPERTY_H
