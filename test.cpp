#include "test.h"
#include <iostream>

void count_atoms(const std::vector<uint8_t>& species) {
    int cnt_Fe = 0, cnt_Gd = 0;
    const int N = static_cast<int>(species.size());

    std::cout << "N = " << N << "\n";
    for (int i = 0; i < N; ++i) {
        if (species[i] == 0) {
            ++cnt_Fe;
        }
        else if (species[i] == 1) {
            ++cnt_Gd;
        }
    }
    std::cout << "Fe count = " << cnt_Fe << "\n";
    std::cout << "Gd count = " << cnt_Gd << "\n";
}