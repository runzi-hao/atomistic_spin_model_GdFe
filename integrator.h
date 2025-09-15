#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <cstdint>
#include <vector>
#include "params.h"

void advance_and_normalize_m(
    const std::vector<double>& mx_in,
    const std::vector<double>& my_in,
    const std::vector<double>& mz_in,
    const std::vector<double>& mx_out,
    const std::vector<double>& my_out,
    const std::vector<double>& mz_out,
    std::vector<double>& dmx_dt,
    std::vector<double>& dmy_dt,
    std::vector<double>& dmz_dt,
    double h_sec);

void advance_and_normalize_m_Heun(
    const std::vector<double>& mx,
    const std::vector<double>& my,
    const std::vector<double>& mz,
    std::vector<double>& dmx_dt_st1,
    std::vector<double>& dmy_dt_st1,
    std::vector<double>& dmz_dt_st1,
    std::vector<double>& dmx_dt_st2,
    std::vector<double>& dmy_dt_st2,
    std::vector<double>& dmz_dt_st2,
    double h_sec);

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
