#ifndef LATTICE_H
#define LATTICE_H
#include <vector>
#include <array>
#include <cstdint>
#include "params.h"

inline int count_fcc_sites(int nx, int ny, int nz) {
    return 4 * nx * ny * nz;
}

/** Build FCC lattice and nearest-neighbor table with PBC. */
void build_fcc_nn(int nx, int ny, int nz,
    std::vector<std::array<int,constants::FCC_NN_COUNT>>& nearest_neighbors);

/** Assign species by fraction (e.g., frac1=0.25 means 25% of sites = 1). */
void assign_species_by_fraction(int N, double frac1,
    std::vector<uint8_t>& species, uint32_t shuffle_seed=0);

#endif //LATTICE_H
