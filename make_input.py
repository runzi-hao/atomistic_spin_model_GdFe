#!/usr/bin/env python3
import argparse, csv, json, random 
from itertools import product
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

# Constants --------------------------------------------------------------------
MU_B = 9.2740100657e-24 # [A m^2]

# Columns
HEADER = [
    "seed", "pre_steps", "run_steps", "save_steps", "dt_sec", 
    "run_parent_path", "run_base_folder", "pre_Te_kelvin", "input_Te_path",
    "nx", "ny", "nz", "a_m", "frac_Gd",
    "J_FeFe_joule_per_link", "J_FeGd_joule_per_link", "J_GdGd_joule_per_link",
    "mx_init_Fe", "my_init_Fe", "mz_init_Fe",
    "mx_init_Gd", "my_init_Gd", "mz_init_Gd",
    "Hx_appl_tesla", "Hy_appl_tesla", "Hz_appl_tesla",
    "mu_ampere_m2_Fe", "alpha_Fe", "gamma_rad_per_tesla_sec_Fe", "ku_joule_per_atom_Fe",
    "easy_axis_x_Fe", "easy_axis_y_Fe", "easy_axis_z_Fe",
    "mu_ampere_m2_Gd", "alpha_Gd", "gamma_rad_per_tesla_sec_Gd", "ku_joule_per_atom_Gd",
    "easy_axis_x_Gd", "easy_axis_y_Gd", "easy_axis_z_Gd"
]

def validate_keys(grid: Dict[str, Any]) -> None:
    missing = [k for k in HEADER if k not in grid]
    extra   = [k for k in grid.keys() if k not in HEADER]
    if missing:
        raise KeyError(f"Missing keys: {missing}")
    if extra:
        raise KeyError(f"Unknown keys: {extra}")

def count_combs(grid: Dict[str, Any]) -> int:
    sizes = [len(grid[k]) for k in HEADER]
    prod = 1
    for s in sizes:
        prod *= s
    return prod

def cartesian_rows(grid: Dict[str, Any]) -> Iterable[List[Any]]:
    value_lists = [grid[k] for k in HEADER]
    for tup in product(*value_lists):
        yield list(tup)

def row_to_map(row: List[Any]) -> Dict[str, Any]:
    return {k: v for k, v in zip(HEADER, row)}

def write_single_csv(path: Path, row: List[Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        w.writerow(row)

def make_timestamp():
    pass

def make_uuid():
    pass

# TODO: Vector.x, .y, .z should not take cartesian product.
# TODO: run_base_folder should be appended by unique identifier.
def main():
    ap = argparse.ArgumentParser(description="Generate input.csv files.")
    ap.add_argument("--config-filename", type=str, 
                    help="JSON file with {key: list}")
    ap.add_argument("--input-filename", type=str, default="input.csv",
                    help="Input filename (default: %(default)s)")
    ap.add_argument("--max-files", type=int, default=100000,
                    help="Restriction on max number of files")
    args = ap.parse_args()

    # Build grid
    if args.config_filename:
        with open(args.config_filename, "r", encoding="utf-8") as f:
            grid = json.load(f)
    else:
        grid = {
            "seed":            [12345],
            "pre_steps":       [40000],
            "run_steps":       [120000],
            "save_steps":      [100],
            "dt_sec":          [0.1e-15],
            "run_parent_path": [r"C:/my_files/research/project_3/project_3_data/data_atomistic_spin_model_GdFe"],
            "run_base_folder": [r"run_20250915_004000_000000"], 
            "pre_Te_kelvin":   [83],
            "input_Te_path":   [r"C:/my_files/research/project_3/project_3_data/data_two_temperature_model/chen_2001/temperatures_20250914_212817_938"],
            "nx":         [50], 
            "ny":         [50], 
            "nz":         [50],
            "a_m":        [1],
            "frac_Gd":    [0.25],
            "J_FeFe_joule_per_link": [2.835e-21],
            "J_FeGd_joule_per_link": [-1.09e-22],
            "J_GdGd_joule_per_link": [1.26e-21],
            "mx_init_Fe":    [0.0], "my_init_Fe":    [0.0], "mz_init_Fe":    [-1.0],
            "mx_init_Gd":    [0.0], "my_init_Gd":    [0.0], "mz_init_Gd":    [1.0],
            "Hx_appl_tesla": [0.0], "Hy_appl_tesla": [0.0], "Hz_appl_tesla": [0.0],
            "mu_ampere_m2_Fe":            [1.92*MU_B], 
            "alpha_Fe":                   [0.05], 
            "gamma_rad_per_tesla_sec_Fe": [1.76e11], 
            "ku_joule_per_atom_Fe":       [0.807246e-23],
            "easy_axis_x_Fe": [0.0], "easy_axis_y_Fe": [0.0], "easy_axis_z_Fe": [1.0],
            "mu_ampere_m2_Gd":            [7.63*MU_B], 
            "alpha_Gd":                   [0.05], 
            "gamma_rad_per_tesla_sec_Gd": [1.76e11], 
            "ku_joule_per_atom_Gd":       [0.807246e-23],
            "easy_axis_x_Gd": [0.0], "easy_axis_y_Gd": [0.0], "easy_axis_z_Gd": [1.0]
        }
    validate_keys(grid)

    total_combs = count_combs(grid)
    if total_combs > args.max_files:
        raise RuntimeError(
            f"Would create --max-files {total_combs} files. "
            "You should reduce list sizes or increase --max-files."
        )
    rows = cartesian_rows(grid) # Build generator

    created  = 0
    for idx, row in enumerate(rows, start=1):
        row_map  = row_to_map(row)
        folder   = Path(row_map["run_parent_path"]) / row_map["run_base_folder"]
        csv_path = folder / args.input_filename
        if csv_path.exists():
            raise FileExistsError(f"File exists: {csv_path}")
        write_single_csv(csv_path, row)
        created += 1
    print(f"Created {created} case(s).")

if __name__ == "__main__":
    main()