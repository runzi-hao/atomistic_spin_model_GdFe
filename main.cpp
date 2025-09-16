#include "params.h"
#include "fields.h"
#include "init.h"
#include "integrator.h"
#include "io.h"
#include "io_temperature_csv.h"
#include "lattice.h"
#include "reductions.h"
#include "rng.h"
#include <filesystem>
#include <iostream>
namespace fs = std::filesystem;

int main() {
    // 1. Read input parameters ------------------------------------------------
    ControlParams control{};
    LatParams lat{};
    MatParams mat[2]{};

    fs::path input_filepath = fs::current_path() / "input.csv";
    read_input_csv(input_filepath.string(), control, lat, mat);
    /// Compute N, normalize easy axes and initial magnetizations
    process_input(lat, mat);

    // 2. Allocate lattice arrays ----------------------------------------------
    // std::vector<double> x_m, y_m, z_m; // site positions
    std::vector<std::array<int, constants::FCC_NN_COUNT>> nearest_neighbors;
    nearest_neighbors.reserve(lat.N);
    std::vector<uint8_t> species; // 0=Fe,1=Gd
    species.reserve(lat.N);

    // 3. Build lattice, assign species ----------------------------------------
    build_fcc_nn(lat.nx, lat.ny, lat.nz, nearest_neighbors);
    RNG rng(control.seed);
    assign_species_by_fraction(lat.N, lat.frac_Gd, species, control.seed);

    // 4. Allocate & initialize other arrays -----------------------------------
    std::vector<double> mx, my, mz;
    std::vector<double> mx_mid(lat.N), my_mid(lat.N), mz_mid(lat.N);
    std::vector<double> dmx_dt_st1(lat.N), dmy_dt_st1(lat.N), dmz_dt_st1(lat.N);
    std::vector<double> dmx_dt_st2(lat.N), dmy_dt_st2(lat.N), dmz_dt_st2(lat.N);
    std::vector<double> Hx_exch_tesla(lat.N), Hy_exch_tesla(lat.N), Hz_exch_tesla(lat.N);
    std::vector<double> Hx_anis_tesla(lat.N), Hy_anis_tesla(lat.N), Hz_anis_tesla(lat.N);
    std::vector<double> Hx_ther_tesla(lat.N), Hy_ther_tesla(lat.N), Hz_ther_tesla(lat.N);
    std::vector<double> Hx_total_tesla(lat.N),Hy_total_tesla(lat.N),Hz_total_tesla(lat.N);

    initialize_m(species, mx, my, mz,
        lat.mx_init_Fe, lat.my_init_Fe, lat.mz_init_Fe,
        lat.mx_init_Gd, lat.my_init_Gd, lat.mz_init_Gd);
    double Hx_appl_tesla=lat.Hx_appl_tesla;
    double Hy_appl_tesla=lat.Hy_appl_tesla;
    double Hz_appl_tesla=lat.Hz_appl_tesla;

    // 5. Read temperature series ----------------------------------------------
    std::vector<double> Te_kelvin_arr;
    read_temperature_series_csv(control.Te_filepath, Te_kelvin_arr);
    if (static_cast<int>(Te_kelvin_arr.size()) < control.run_steps) {
        throw std::runtime_error("Temperature data not enough (" +
            std::to_string(control.run_steps) + " required)");
    }

    // 6. Time evolve ----------------------------------------------------------
    for (int curr_step=0; curr_step <= control.pre_steps + control.run_steps;
        ++curr_step)
    {
        // Set temperature
        double T_kelvin = (curr_step < control.pre_steps) ?
            control.pre_Te_kelvin :
            Te_kelvin_arr[curr_step - control.pre_steps];
        compute_ther_field_once(mat, species, T_kelvin, control.dt_sec, rng,
            Hx_ther_tesla, Hy_ther_tesla, Hz_ther_tesla);

        // Heun stage-1 --------------------------------------------------------
        advance_and_normalize_m(mx, my, mz,
            mx_mid, my_mid, mz_mid,
            dmx_dt_st1, dmy_dt_st1, dmz_dt_st1, 0.0);
        compute_total_field(mat, lat.J_joule_per_link, nearest_neighbors,
            species,
            mx_mid, my_mid, mz_mid,
            Hx_appl_tesla,  Hy_appl_tesla,  Hz_appl_tesla,
            Hx_exch_tesla,  Hy_exch_tesla,  Hz_exch_tesla,
            Hx_anis_tesla,  Hy_anis_tesla,  Hz_anis_tesla,
            Hx_ther_tesla,  Hy_ther_tesla,  Hz_ther_tesla,
            Hx_total_tesla, Hy_total_tesla, Hz_total_tesla);
        compute_dm_dt_kernel(mat, species,
            mx_mid, my_mid, mz_mid,
            Hx_total_tesla, Hy_total_tesla, Hz_total_tesla,
            dmx_dt_st1, dmy_dt_st1, dmz_dt_st1);

        // Heun stage-2 --------------------------------------------------------
        advance_and_normalize_m(mx, my, mz,
            mx_mid, my_mid, mz_mid,
            dmx_dt_st1, dmy_dt_st1, dmz_dt_st1, control.dt_sec);
        compute_total_field(mat, lat.J_joule_per_link, nearest_neighbors,
            species,
            mx_mid, my_mid, mz_mid,
            Hx_appl_tesla,  Hy_appl_tesla,  Hz_appl_tesla,
            Hx_exch_tesla,  Hy_exch_tesla,  Hz_exch_tesla,
            Hx_anis_tesla,  Hy_anis_tesla,  Hz_anis_tesla,
            Hx_ther_tesla,  Hy_ther_tesla,  Hz_ther_tesla,
            Hx_total_tesla, Hy_total_tesla, Hz_total_tesla);
        compute_dm_dt_kernel(mat, species,
            mx_mid, my_mid, mz_mid,
            Hx_total_tesla, Hy_total_tesla, Hz_total_tesla,
            dmx_dt_st2, dmy_dt_st2, dmz_dt_st2);

        // Advance m -----------------------------------------------------------
        advance_and_normalize_m_Heun(mx, my, mz,
            dmx_dt_st1, dmy_dt_st1, dmz_dt_st1,
            dmx_dt_st2, dmy_dt_st2, dmz_dt_st2, control.dt_sec);

        // 7. Reductions & outputs ---------------------------------------------
        if (curr_step % control.save_steps == 0) {
            BulkValues bulk_vals{};
            compute_bulk_m(species, mx, my, mz, bulk_vals);
            fs::path run_dir = fs::path(control.run_parent_dir) /
                control.run_base_folder;
            fs::path run_filepath = run_dir / "bulk_values_vs_time.csv";
            write_bulk_values(run_filepath.string(), curr_step, T_kelvin, bulk_vals);
        }
        if (curr_step % control.show_steps == 0) {
            std::cout << curr_step << std::endl;
        }
    }

    return 0;
}
