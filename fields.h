#ifndef FIELDS_H
#define FIELDS_H
#include <cstdint>
#include <vector>
#include "params.h"
#include "rng.h"

void compute_exch_field(
    const double J_joule_per_link[2][2],
    const std::vector<std::array<int, constants::FCC_NN_COUNT>>&
    nearest_neighbors,
    const std::vector<uint8_t>& species,
    const std::vector<double>& mx,
    const std::vector<double>& my,
    const std::vector<double>& mz,
    std::vector<double>& Hx_exch_tesla,
    std::vector<double>& Hy_exch_tesla,
    std::vector<double>& Hz_exch_tesla);

void compute_uniaxial_anis_field(
    const MatParams phys_params[2],
    const std::vector<uint8_t>& species,
    const std::vector<double>& mx_arr,
    const std::vector<double>& my_arr,
    const std::vector<double>& mz_arr,
    std::vector<double>& Hx_anis_tesla,
    std::vector<double>& Hy_anis_tesla,
    std::vector<double>& Hz_anis_tesla);

void compute_ther_field_once(
    const MatParams mat[2],
    const std::vector<uint8_t>& species,
    double T_kelvin, double dt_sec, RNG& rng,
    std::vector<double>& Hx_ther_tesla,
    std::vector<double>& Hy_ther_tesla,
    std::vector<double>& Hz_ther_tesla);

void compute_total_field(
    const MatParams phys_params[2],
    const double J_joule_per_link[2][2],
    const std::vector<std::array<int, constants::FCC_NN_COUNT>>&
    nearest_neighbors,
    const std::vector<uint8_t>& species,
    const std::vector<double>& mx,
    const std::vector<double>& my,
    const std::vector<double>& mz,
    double Hx_appl_tesla, double Hy_appl_tesla, double Hz_appl_tesla,
    std::vector<double>& Hx_exch_tesla,
    std::vector<double>& Hy_exch_tesla,
    std::vector<double>& Hz_exch_tesla,
    std::vector<double>& Hx_anis_tesla,
    std::vector<double>& Hy_anis_tesla,
    std::vector<double>& Hz_anis_tesla,
    const std::vector<double>& Hx_ther_tesla,
    const std::vector<double>& Hy_ther_tesla,
    const std::vector<double>& Hz_ther_tesla,
    std::vector<double>& Hx_total_tesla,
    std::vector<double>& Hy_total_tesla,
    std::vector<double>& Hz_total_tesla);

#endif //FIELDS_H
