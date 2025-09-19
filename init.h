#ifndef INIT_H
#define INIT_H
#include <cassert>
#include <cstdint>
#include <vector>

/// Initialize per-site magnetization:
inline void initialize_m(const uint8_t* species, int N,
    double* mx, double* my, double* mz,
    const double mx_Fe, const double my_Fe, const double mz_Fe,
    const double mx_Gd, const double my_Gd, const double mz_Gd)
{
    assert(species && mx && my && mz);
    for (int i = 0; i < N; ++i) {
        const bool isFe = (species[i] == 0);
        mx[i] = isFe ? mx_Fe : mx_Gd;
        my[i] = isFe ? my_Fe : my_Gd;
        mz[i] = isFe ? mz_Fe : mz_Gd;
    }
}

/// Convenience wrapper for std::vector storage.
/// Resizes mx/my/mz to match species.size() if needed.
inline void initialize_m(const std::vector<uint8_t>& species,
    std::vector<double>& mx,
    std::vector<double>& my,
    std::vector<double>& mz,
    double mx_Fe, double my_Fe, double mz_Fe,
    double mx_Gd, double my_Gd, double mz_Gd)
{
    const int N = static_cast<int>(species.size());
    mx.resize(N); my.resize(N); mz.resize(N);
    initialize_m(species.data(), N,
        mx.data(), my.data(), mz.data(),
        mx_Fe, my_Fe,  mz_Fe,
        mx_Gd, my_Gd, mz_Gd);
}

#endif //INIT_H
