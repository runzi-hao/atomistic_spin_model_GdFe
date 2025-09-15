#include "io.h"
#include "io_csv_utils.h"
#include "lattice.h"
#include "math_utils.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace io::csv;

/// Expected columns (order free):
// run_parent_dir, run_base_folder, Te_filepath,
// seed, pre_steps, run_steps, save_steps, show_steps, dt_sec, pre_Te_kelvin,
// nx, ny, nz, a_m, frac_Gd,
// J_FeFe_joule_per_link, J_FeGd_joule_per_link, J_GdGd_joule_per_link,
// mx_init_Fe, my_init_Fe, mz_init_Fe,
// mx_init_Gd, my_init_Gd, mz_init_Gd,
// Hx_appl_tesla, Hy_appl_tesla, Hz_appl_tesla,
// mu_ampere_m2_Fe, alpha_Fe, gamma_rad_per_tesla_sec_Fe, ku_joule_per_atom_Fe,
// easy_axis_x_Fe, easy_axis_y_Fe, easy_axis_z_Fe
// mu_ampere_m2_Gd, alpha_Gd, gamma_rad_per_tesla_sec_Gd, ku_joule_per_atom_Gd,
// easy_axis_x_Gd, easy_axis_y_Gd, easy_axis_z_Gd

namespace {
    double get_dou(const std::unordered_map<std::string, int>& idx,
        const std::vector<std::string>& vals, const std::string& key)
    {
        auto it = idx.find(key);
        if (it == idx.end())
            throw std::runtime_error(std::string("Missing key: ") + key);
        return std::stod(vals[it->second]);
    }
    int get_int(const std::unordered_map<std::string,int>& idx,
        const std::vector<std::string>& vals, const std::string& key)
    {
        auto it = idx.find(key);
        if (it == idx.end())
            throw std::runtime_error(std::string("Missing key: ") + key);
        return std::stoi(vals[it->second]);
    }
    uint32_t get_u32(const std::unordered_map<std::string,int>& idx,
        const std::vector<std::string>& vals, const std::string& key)
    {
        auto it = idx.find(key);
        if (it == idx.end())
            throw std::runtime_error(std::string("Missing key: ") + key);
        return static_cast<uint32_t>(std::stoul(vals[it->second]));
    }
    std::string get_str(const std::unordered_map<std::string,int>& idx,
        const std::vector<std::string>& vals, const std::string& key)
    {
        auto it = idx.find(key);
        if (it == idx.end())
            throw std::runtime_error(std::string("Missing key: ") + key);
        return vals[it->second];
    }
}

bool read_input_csv(const std::string& csv_path, ControlParams& control,
    LatParams& lat, MatParams mat[2])
{
    std::ifstream fin(csv_path);
    if (!fin) {
        std::cerr << "io: Failed to open file: " << csv_path << "\n";
        return false;
    }

    std::string header, values;
    if (!std::getline(fin, header) || !std::getline(fin, values)) {
        std::cerr << "io: Expected 2 lines (header + values).\n";
        return false;
    }

    auto keys = split_line_csv(header);
    auto vals_str = split_line_csv(values);
    if (vals_str.size() != keys.size()) {
        std::cerr << "io: Expected " << keys.size() << " entries for values.\n";
        return false;
    }

    std::unordered_map<std::string,int> key_idx_map;
    for (int i=0; i < static_cast<int>(keys.size()); ++i) {
        if (const std::string& k = keys[i]; !k.empty())
            key_idx_map[k] = i;
    }

    try {
        // Sim params
        control.seed            = get_u32(key_idx_map, vals_str, "seed");
        control.pre_steps       = get_int(key_idx_map, vals_str, "pre_steps");
        control.run_steps       = get_int(key_idx_map, vals_str, "run_steps");
        control.save_steps      = get_int(key_idx_map, vals_str, "save_steps");
        control.show_steps      = get_int(key_idx_map, vals_str, "show_steps");
        control.dt_sec          = get_dou(key_idx_map, vals_str, "dt_sec");
        control.pre_Te_kelvin   = get_dou(key_idx_map, vals_str, "pre_Te_kelvin");
        control.run_parent_dir = get_str(key_idx_map, vals_str, "run_parent_dir");
        control.run_base_folder = get_str(key_idx_map, vals_str, "run_base_folder");
        control.Te_filepath   = get_str(key_idx_map, vals_str, "Te_filepath");

        // FCC params
        lat.nx      = get_int(key_idx_map, vals_str, "nx");
        lat.ny      = get_int(key_idx_map, vals_str, "ny");
        lat.nz      = get_int(key_idx_map, vals_str, "nz");
        lat.a_m     = get_dou(key_idx_map, vals_str, "a_m");
        lat.frac_Gd = get_dou(key_idx_map, vals_str, "frac_Gd");
        {
            const double J_FeFe_joule_per_link = get_dou(key_idx_map, vals_str,
                "J_FeFe_joule_per_link");
            const double J_FeGd_joule_per_link = get_dou(key_idx_map, vals_str,
                "J_FeGd_joule_per_link");
            const double J_GdGd_joule_per_link = get_dou(key_idx_map, vals_str,
                "J_GdGd_joule_per_link");
            lat.J_joule_per_link[0][0] = J_FeFe_joule_per_link;
            lat.J_joule_per_link[0][1] = J_FeGd_joule_per_link;
            lat.J_joule_per_link[1][0] = J_FeGd_joule_per_link;
            lat.J_joule_per_link[1][1] = J_GdGd_joule_per_link;
        }
        lat.mx_init_Fe  = get_dou(key_idx_map, vals_str, "mx_init_Fe");
        lat.my_init_Fe  = get_dou(key_idx_map, vals_str, "my_init_Fe");
        lat.mz_init_Fe  = get_dou(key_idx_map, vals_str, "mz_init_Fe");
        lat.mx_init_Gd  = get_dou(key_idx_map, vals_str, "mx_init_Gd");
        lat.my_init_Gd  = get_dou(key_idx_map, vals_str, "my_init_Gd");
        lat.mz_init_Gd  = get_dou(key_idx_map, vals_str, "mz_init_Gd");
        lat.Hx_appl_tesla = get_dou(key_idx_map, vals_str, "Hx_appl_tesla");
        lat.Hy_appl_tesla = get_dou(key_idx_map, vals_str, "Hy_appl_tesla");
        lat.Hz_appl_tesla = get_dou(key_idx_map, vals_str, "Hz_appl_tesla");

        // Phys params by species
        auto fill_phys = [&mat, &key_idx_map, &vals_str](int s, const char* tag)
        {
            MatParams& P = mat[s];
            const std::string mu = std::string("mu_ampere_m2_") + tag;
            const std::string al = std::string("alpha_") + tag;
            const std::string ga = std::string("gamma_rad_per_tesla_sec_") + tag;
            const std::string ku = std::string("ku_joule_per_atom_")    + tag;
            const std::string ea_x = std::string("easy_axis_x_") + tag;
            const std::string ea_y = std::string("easy_axis_y_") + tag;
            const std::string ea_z = std::string("easy_axis_z_") + tag;

            P.mu_ampere_m2            = get_dou(key_idx_map, vals_str, mu);
            P.alpha                   = get_dou(key_idx_map, vals_str, al);
            P.gamma_rad_per_tesla_sec = get_dou(key_idx_map, vals_str, ga);
            P.ku_joule_per_atom       = get_dou(key_idx_map, vals_str, ku);
            P.easy_axis.x = get_dou(key_idx_map, vals_str, ea_x);
            P.easy_axis.y = get_dou(key_idx_map, vals_str, ea_y);
            P.easy_axis.z = get_dou(key_idx_map, vals_str, ea_z);
        };
        fill_phys(0, "Fe");
        fill_phys(1, "Gd");
    }
    catch (const std::exception& e) {
        std::cerr << "io: " << e.what() << "\n";
        return false;
    }
    return true;
}

void process_input(LatParams& lat, MatParams mat[2]) {
    lat.N = count_fcc_sites(lat.nx, lat.ny, lat.nz);
    normalize3(lat.mx_init_Fe, lat.my_init_Fe, lat.mz_init_Fe);
    normalize3(lat.mx_init_Gd, lat.my_init_Gd, lat.mz_init_Gd);
    for (int i=0; i < 2; ++i) {
        auto& m = mat[i];
        normalize3(m.easy_axis);
    }
}

void write_bulk_values(const std::string& csv_path, const int time_step,
    const double T_kelvin, const BulkValues& bulk)
{
    try {
        std::filesystem::path p(csv_path);
        if (!p.parent_path().empty()) {
            std::filesystem::create_directories(p.parent_path());
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string(
            "write_bulk_values: Failed to create directories: ") + e.what());
    }

    // Check if header is needed
    bool need_header = false;
    if (!std::filesystem::exists(csv_path) ||
        std::filesystem::file_size(csv_path) == 0)
    {
        need_header = true;
    }

    std::ofstream ofs(csv_path, std::ios::out | std::ios::app);
    if (!ofs) {
        throw std::runtime_error(
            "write_bulk_values: Failed to open file: " + csv_path);
    }
    ofs.imbue(std::locale::classic());
    ofs << std::defaultfloat << std::setprecision(10);

    if (need_header) {
        ofs << "time_step,T_kelvin,"
               "mx_Fe,my_Fe,mz_Fe,"
               "mx_Gd,my_Gd,mz_Gd,"
               "mx_bulk,my_bulk,mz_bulk\n";
    }

    ofs
        << time_step << ',' << T_kelvin << ','
        << bulk.mx_Fe   << ',' << bulk.my_Fe   << ',' << bulk.mz_Fe << ','
        << bulk.mx_Gd   << ',' << bulk.my_Gd   << ',' << bulk.mz_Gd << ','
        << bulk.mx_bulk << ',' << bulk.my_bulk << ',' << bulk.mz_bulk
        << '\n';
}


