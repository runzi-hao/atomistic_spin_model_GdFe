#ifndef IO_H
#define IO_H
#include <string>
#include "params.h"

bool read_input_csv(const std::string& csv_path, ControlParams& control,
    LatParams& lat, MatParams mat[2]);

void process_input(LatParams& lat, MatParams mat[2]);

void write_bulk_values(const std::string& csv_path, int time_step,
    double T_kelvin, const BulkValues& bulk);

#endif //IO_H
