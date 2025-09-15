#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <cstdint>
#include <vector>
#include "params.h"

void advance_and_normalize_m(
    const std::vector<double>& mx_in,
    const std::vector<double>& my_in,
    const std::vector<double>& mz_in,
    std::vector<double>& mx_out,
    std::vector<double>& my_out,
    std::vector<double>& mz_out,
    const std::vector<double>& dmx_dt,
    const std::vector<double>& dmy_dt,
    const std::vector<double>& dmz_dt,
    double h_sec);

void advance_and_normalize_m_Heun(
    std::vector<double>& mx,
    std::vector<double>& my,
    std::vector<double>& mz,
    const std::vector<double>& dmx_dt_st1,
    const std::vector<double>& dmy_dt_st1,
    const std::vector<double>& dmz_dt_st1,
    const std::vector<double>& dmx_dt_st2,
    const std::vector<double>& dmy_dt_st2,
    const std::vector<double>& dmz_dt_st2,
    const double h_sec);

void compute_dm_dt_kernel(
    const MatParams phys_params[2],
    const std::vector<uint8_t>& species,
    const std::vector<double>& mx_arr,
    const std::vector<double>& my_arr,
    const std::vector<double>& mz_arr,
    const std::vector<double>& Hx_total_tesla_arr,
    const std::vector<double>& Hy_total_tesla_arr,
    const std::vector<double>& Hz_total_tesla_arr,
    std::vector<double>& dmx_dt,
    std::vector<double>& dmy_dt,
    std::vector<double>& dmz_dt);

#endif //INTEGRATOR_H
