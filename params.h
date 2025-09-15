#ifndef PARAMS_H
#define PARAMS_H
#include <cstdint>
#include <string>

struct Vec3 {
    double x, y, z;
};
namespace constants {
    constexpr double KB_JOULE_PER_KELVIN = 1.380649e-23;
    constexpr int FCC_NN_COUNT    = 12;
    constexpr int FCC_BASIS_COUNT = 4;
}

//------------------------------------------------------------------------------
struct ControlParams {
    uint32_t seed;
    int pre_steps;
    int run_steps;
    int save_steps;
    double dt_sec;
    double pre_Te_kelvin;
    std::string run_parent_path;
    std::string run_base_folder;
    std::string input_Te_path;
};
struct LatParams {
    int nx, ny, nz; // number of cells
    double a_m;
    double frac_Gd;
    double J_joule_per_link[2][2]; // {{J_FeFe, J_FeGd}, {J_FeGd, J_GdGd}}
    double mx_init_Fe, my_init_Fe, mz_init_Fe;
    double mx_init_Gd, my_init_Gd, mz_init_Gd;
    double Hx_appl_tesla, Hy_appl_tesla, Hz_appl_tesla;
    /// Derived
    int N;
};
struct MatParams {
    double mu_ampere_m2;
    double alpha;
    double gamma_rad_per_tesla_sec;
    double ku_joule_per_atom;
    Vec3 easy_axis;
};

//------------------------------------------------------------------------------
struct BulkValues {
    double mx_Fe,   my_Fe,   mz_Fe;
    double mx_Gd,   my_Gd,   mz_Gd;
    double mx_bulk, my_bulk, mz_bulk;
};

// // std::vector<double> x_m, y_m, z_m; // site positions
// std::vector<std::array<int, constants::FCC_NN_COUNT>> nearest_neighbors;
// std::vector<uint8_t> species;         // 0=Fe,1=Gd
// // Magnetic parameters
// std::vector<double> mx, my, mz;
// std::vector<double> mx_mid, my_mid, mz_mid;
// std::vector<double> dmx_dt_st1, dmy_dt_st1, dmz_dt_st1;
// std::vector<double> dmx_dt_st2, dmy_dt_st2, dmz_dt_st2;
// std::vector<double> Hx_exch_tesla,  Hy_exch_tesla,  Hz_exch_tesla;
// std::vector<double> Hx_anis_tesla,  Hy_anis_tesla,  Hz_anis_tesla;
// std::vector<double> Hx_ther_tesla,  Hy_ther_tesla,  Hz_ther_tesla;
// std::vector<double> Hx_total_tesla, Hy_total_tesla, Hz_total_tesla;

#endif //PARAMS_H
