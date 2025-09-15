#include "integrator.h"
#include "math_utils.h"
#include <vector>

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
    const double h_sec)
{
    for (int i=0; i < mx_in.size(); ++i) {
        mx_out[i] = mx_in[i] + h_sec * dmx_dt[i];
        my_out[i] = my_in[i] + h_sec * dmy_dt[i];
        mz_out[i] = mz_in[i] + h_sec * dmz_dt[i];
    }
    for (int i=0; i < mx_out.size(); ++i) {
        normalize3(mx_out[i], my_out[i], mz_out[i]);
    }
}

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
    const double h_sec)
{
    for (int i=0; i < mx.size(); ++i) {
        mx[i] += h_sec * 0.5 * (dmx_dt_st1[i] + dmx_dt_st2[i]);
        my[i] += h_sec * 0.5 * (dmy_dt_st1[i] + dmy_dt_st2[i]);
        mz[i] += h_sec * 0.5 * (dmz_dt_st1[i] + dmz_dt_st2[i]);
        normalize3(mx[i], my[i], mz[i]);
    }
}

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
    std::vector<double>& dmz_dt)
{
    for (int i=0; i < species.size(); ++i) {
        const int s = species[i];
        const double gamma_rad_per_tesla_sec =
            phys_params[s].gamma_rad_per_tesla_sec;
        const double alpha                   = phys_params[s].alpha;
        const double mx = mx_arr[i];
        const double my = my_arr[i];
        const double mz = mz_arr[i];
        const double Hx_total_tesla = Hx_total_tesla_arr[i];
        const double Hy_total_tesla = Hy_total_tesla_arr[i];
        const double Hz_total_tesla = Hz_total_tesla_arr[i];

        const double gamma_prime_rad_per_tesla_sec =
            -gamma_rad_per_tesla_sec / (1. + alpha*alpha);

        // c1 = m x H
        const double c1x = my*Hz_total_tesla - mz*Hy_total_tesla;
        const double c1y = mz*Hx_total_tesla - mx*Hz_total_tesla;
        const double c1z = mx*Hy_total_tesla - my*Hx_total_tesla;

        // c2 = m x (m x H)
        const double c2x = my*c1z - mz*c1y;
        const double c2y = mz*c1x - mx*c1z;
        const double c2z = mx*c1y - my*c1x;

        dmx_dt[i] = gamma_prime_rad_per_tesla_sec * ( c1x + alpha * c2x );
        dmy_dt[i] = gamma_prime_rad_per_tesla_sec * ( c1y + alpha * c2y );
        dmz_dt[i] = gamma_prime_rad_per_tesla_sec * ( c1z + alpha * c2z );
    }
}