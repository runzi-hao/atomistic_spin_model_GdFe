#ifndef IO_TEMPERATURE_CSV_H
#define IO_TEMPERATURE_CSV_H
#include <string>
#include "temperature_series.h"

/** Read CSV with only one column Te_kelvin */
void read_temperature_series_csv(const std::string& path,
    std::vector<double>& Te_kelvin);

#endif //IO_TEMPERATURE_CSV_H
