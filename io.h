#ifndef IO_H
#define IO_H
#include <string>
#include <vector>
#include "params.h"

bool read_input_csv(const std::string& csv_path, ControlParams& control,
    LatParams& lat, MatParams mat[2]);

void process_input(LatParams& lat, MatParams mat[2]);

void write_bulk_values(const std::string& csv_path, int time_step,
    double T_kelvin, const BulkValues& bulk);

void write_bulk_values(const std::string& csv_path, int time_step,
    double T_kelvin,
    const BulkValues& bulk_vals,
    const BulkFields& bulk_fields);

void write_nearest_neighbors(const std::string& filepath,
    int nx, int ny, int nz, int n_basis,
    const std::vector<std::array<int, constants::FCC_NN_COUNT>>&
    nearest_neighbors);

void write_site_species(const std::string& filepath,
    int nx, int ny, int nz, int n_basis,
    const std::vector<uint8_t>& species);

#endif //IO_H
