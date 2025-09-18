#ifndef LATTICE_H
#define LATTICE_H
#include <vector>
#include <array>
#include <cstdint>
#include "params.h"

inline int count_fcc_sites(const int nx, const int ny, const int nz) {
    return 4 * nx * ny * nz;
}

// Invert p = (((i*ny + j)*nz + k)*n_basis + b)
inline void invert_linear_index(int p,
    const int nx, const int ny, const int nz, const int n_basis,
    int& i, int& j, int& k, int& b)
{
    b = p % n_basis;    p /= n_basis;
    k = p % nz;         p /= nz;
    j = p % ny;         p /= ny;
    i = p;
}

/** Build FCC lattice and nearest-neighbor table with PBC. */
void build_fcc_nn(int nx, int ny, int nz,
    std::vector<std::array<int,constants::FCC_NN_COUNT>>& nearest_neighbors);

/** Assign species by fraction (e.g., frac1=0.25 means 25% of sites = 1). */
void assign_species_by_fraction(int N, double frac1,
    std::vector<uint8_t>& species, uint32_t shuffle_seed=0);

#endif //LATTICE_H
