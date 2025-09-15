#include "lattice.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <random>

namespace {
    // Basis offsets in half-steps
    constexpr int BASIS_OFF[4][3] = {
        {0,0,0}, {0,1,1}, {1,0,1}, {1,1,0}
    };

    // Neighbor offsets (in half-steps):
    // 12 vectors where two axes are ±1 and one is 0.
    constexpr int NN12_OFF[12][3] = {
        { 0,  1,  1}, { 0,  1, -1}, { 0, -1,  1}, { 0, -1, -1},
        { 1,  0,  1}, { 1,  0, -1}, {-1,  0,  1}, {-1,  0, -1},
        { 1,  1,  0}, { 1, -1,  0}, {-1,  1,  0}, {-1, -1,  0}
    };

    // FCC basis offsets (in half-steps):
    // b=0:(0,0,0), b=1:(0,1,1), b=2:(1,0,1), b=3:(1,1,0).
    int fcc_basis_from_remainders(const int rx, const int ry, const int rz) {
        if (rx==0 && ry==0 && rz==0) return 0;
        if (rx==0 && ry==1 && rz==1) return 1;
        if (rx==1 && ry==0 && rz==1) return 2;
        if (rx==1 && ry==1 && rz==0) return 3;
        assert(false && "Invalid (rx,ry,rz) for FCC basis, expected from "
                        "{000, 011, 101, 110})");
        return -1;
    }

    // Linear index from (cell i,j,k) and basis b
    int lin_from_cell_and_basis(
        const int i, const int j, const int k, const int b,
        const int nx, const int ny, const int nz, const int n_basis)
    {
        return (((i*ny + j)*nz + k)*n_basis + b);
    }

    //  10 %  3 =  1
    //  10 % -3 =  1
    // -10 %  3 = -1
    // -10 % -3 = -1
    int wrap_modulo(int v, int m) {
        assert(m > 0 && "wrap_modulo: modulus m must be > 0");
        int r = v % m;
        return (r < 0) ? (r + m) : r; // Returns in [0, m-1]
    }
}

void build_fcc_nn(const int nx, const int ny, const int nz,
    std::vector<std::array<int,constants::FCC_NN_COUNT>>& nearest_neighbors)
{
    const int N  = count_fcc_sites(nx, ny, nz);
    const int GX = 2*nx, GY = 2*ny, GZ = 2*nz; // doubled-grid extents

    nearest_neighbors.resize(N);

    // For each site (i,j,k,b) with grid coords (gx,gy,gz) =
    // (2i+bx, 2j+by, 2k+bz), add the 12 NN integer offsets, wrap with PBC, map
    // back to (cell, basis).
    int p = 0;
    for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            for (int k=0; k<nz; ++k) {
                for (int b=0; b < constants::FCC_BASIS_COUNT; ++b, ++p) {
                    const int gx = 2*i + BASIS_OFF[b][0];
                    const int gy = 2*j + BASIS_OFF[b][1];
                    const int gz = 2*k + BASIS_OFF[b][2];

                    std::array<int,constants::FCC_NN_COUNT> neis;
                    for (int q=0; q < constants::FCC_NN_COUNT; ++q) {
                        int nei_gx = wrap_modulo(gx + NN12_OFF[q][0], GX);
                        int nei_gy = wrap_modulo(gy + NN12_OFF[q][1], GY);
                        int nei_gz = wrap_modulo(gz + NN12_OFF[q][2], GZ);

                        // grid coords back to cell + basis:
                        int nei_i = nei_gx >> 1;  // divide by 2
                        int nei_j = nei_gy >> 1;
                        int nei_k = nei_gz >> 1;
                        int rx = nei_gx & 1;      // remainder (0 or 1)
                        int ry = nei_gy & 1;
                        int rz = nei_gz & 1;
                        int nei_b = fcc_basis_from_remainders(rx, ry, rz);

                        neis[q] = lin_from_cell_and_basis(nei_i, nei_j, nei_k,
                            nei_b, nx, ny, nz, constants::FCC_BASIS_COUNT);
                    }
                    nearest_neighbors[p] = neis;
                }
            }
        }
    }
}

// Reproducible random assignment by fraction using hash+sort trick
void assign_species_by_fraction(const int N, double frac1,
    std::vector<uint8_t>& species, uint32_t shuffle_seed)
{
    species.resize(N, 0);
    if (N == 0) return;

    frac1 = std::clamp(frac1, 0.0, 1.0);
    const int target_ones = static_cast<int>(std::round(frac1 * N));
    if (target_ones <= 0) return;
    if (target_ones == N) {species.assign(N, 1); return;}

    std::vector<uint64_t> keys(N);
    std::vector<int>      idx(N);
    std::iota(idx.begin(), idx.end(), 0);

    // splitmix64 mixer
    auto mix = [shuffle_seed](uint64_t x){
        x += (uint64_t)0x9E3779B97F4A7C15ULL + shuffle_seed;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        x =  x ^ (x >> 31);
        return x;
    };
    for (int i=0; i<N; ++i) keys[i] = mix(static_cast<uint64_t>(i));

    std::sort(idx.begin(), idx.end(),
        [&keys](const int a, const int b){ return keys[a] < keys[b]; });

    for (int r=0; r<target_ones; ++r) species[idx[r]] = 1;
}
