#ifndef TEMPERATURE_SERIES_H
#define TEMPERATURE_SERIES_H
#include <vector>
#include <stdexcept>
#include "math_utils.h"

struct TemperatureSeries {
    std::vector<int> time_step;
    std::vector<double> time_sec, Te_kelvin, Tp_kelvin; // Tp_kelvin optional

    bool   has_lattice_temperature() const { return !Tp_kelvin.empty(); }
    double sample_Te(double t) const { return interpolate_clamped(time_sec, Te_kelvin, t); }
    double sample_Tp(double t) const {
        if (!has_lattice_temperature())
            throw std::runtime_error("Tp_kelvin not present");
        return interpolate_clamped(time_sec, Tp_kelvin, t);
    }
};

#endif //TEMPERATURE_SERIES_H
