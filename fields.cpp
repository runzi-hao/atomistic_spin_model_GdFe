#include "fields.h"
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

void compute_exch_field(
    const MatParams mat[2],
    const double J_joule_per_link[2][2],
    const std::vector<std::array<int, constants::FCC_NN_COUNT>>&
    nearest_neighbors,
    const std::vector<uint8_t>& species,
    const std::vector<double>& mx,
    const std::vector<double>& my,
    const std::vector<double>& mz,
    std::vector<double>& Hx_exch_tesla,
    std::vector<double>& Hy_exch_tesla,
    std::vector<double>& Hz_exch_tesla)
{
    const int N = static_cast<int>(species.size());
    for (int i=0; i < N; ++i) {
        const int si = species[i];
        const double mu_ampere_m2 = mat[si].mu_ampere_m2;
        Hx_exch_tesla[i] = 0.0;
        Hy_exch_tesla[i] = 0.0;
        Hz_exch_tesla[i] = 0.0;
        for (int k=0; k < nearest_neighbors[i].size(); ++k) {
            int j = nearest_neighbors[i][k];
            const int sj = species[j];
            // TODO: Discuss again whether should use factor of 2.
            const double J_ij_joule_per_link = J_joule_per_link[si][sj] * 1.0;;
            Hx_exch_tesla[i] += J_ij_joule_per_link * mx[j] / mu_ampere_m2;
            Hy_exch_tesla[i] += J_ij_joule_per_link * my[j] / mu_ampere_m2;
            Hz_exch_tesla[i] += J_ij_joule_per_link * mz[j] / mu_ampere_m2;
        }
    }
}

void compute_uniaxial_anis_field(
    const MatParams mat[2],
    const std::vector<uint8_t>& species,
    const std::vector<double>& mx_arr,
    const std::vector<double>& my_arr,
    const std::vector<double>& mz_arr,
    std::vector<double>& Hx_anis_tesla,
    std::vector<double>& Hy_anis_tesla,
    std::vector<double>& Hz_anis_tesla)
{
    const int N = static_cast<int>(species.size());
    for (int i=0; i < N; ++i) {
        const int s = species[i];
        const double mu_ampere_m2      = mat[s].mu_ampere_m2;
        const double ku_joule_per_atom = mat[s].ku_joule_per_atom;
        const Vec3   easy_axis         = mat[s].easy_axis;
        const double mx = mx_arr[i];
        const double my = my_arr[i];
        const double mz = mz_arr[i];

        const double dot =
            mx*easy_axis.x + my*easy_axis.y + mz*easy_axis.z;
        Hx_anis_tesla[i] = 2.*ku_joule_per_atom*dot*easy_axis.x / mu_ampere_m2;
        Hy_anis_tesla[i] = 2.*ku_joule_per_atom*dot*easy_axis.y / mu_ampere_m2;
        Hz_anis_tesla[i] = 2.*ku_joule_per_atom*dot*easy_axis.z / mu_ampere_m2;
    }
}

// TODO: May use one normal_distribution instead of building three in each step.
void compute_ther_field_once(const MatParams mat[2],
    const std::vector<uint8_t>& species,
    double T_kelvin, double dt_sec, RNG& rng,
    std::vector<double>& Hx_ther_tesla,
    std::vector<double>& Hy_ther_tesla,
    std::vector<double>& Hz_ther_tesla)
{
    const int N = static_cast<int>(species.size());
    for (int i=0; i < N; ++i) {
        const int s = species[i];
        const double alpha = mat[s].alpha;
        const double gamma_rad_per_tesla_sec =
            mat[s].gamma_rad_per_tesla_sec;
        const double mu_ampere_m2 = mat[s].mu_ampere_m2;

        const double sigma_tesla =
            std::sqrt(2. * alpha * constants::KB_JOULE_PER_KELVIN * T_kelvin
                / (gamma_rad_per_tesla_sec * mu_ampere_m2 * dt_sec));
        Hx_ther_tesla[i] = sigma_tesla * rng.normal(0.0, 1.0);
        Hy_ther_tesla[i] = sigma_tesla * rng.normal(0.0, 1.0);
        Hz_ther_tesla[i] = sigma_tesla * rng.normal(0.0, 1.0);
    }
}

void compute_total_field(
    const MatParams mat[2],
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
    std::vector<double>& Hz_total_tesla)
{
    compute_exch_field(mat, J_joule_per_link, nearest_neighbors, species,
        mx, my, mz,
        Hx_exch_tesla, Hy_exch_tesla, Hz_exch_tesla);
    compute_uniaxial_anis_field(mat, species,
        mx, my, mz,
        Hx_anis_tesla, Hy_anis_tesla, Hz_anis_tesla);

    for (int i=0; i < Hx_total_tesla.size(); ++i) {
        Hx_total_tesla[i] = Hx_appl_tesla +
            Hx_exch_tesla[i] + Hx_anis_tesla[i] + Hx_ther_tesla[i];
        Hy_total_tesla[i] = Hy_appl_tesla +
            Hy_exch_tesla[i] + Hy_anis_tesla[i] + Hy_ther_tesla[i];
        Hz_total_tesla[i] = Hz_appl_tesla +
            Hz_exch_tesla[i] + Hz_anis_tesla[i] + Hz_ther_tesla[i];
    }
}