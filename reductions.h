#ifndef REDUCTIONS_H
#define REDUCTIONS_H
#include <cstdint>
#include <vector>
#include "params.h"

/// species[i]: 0 = Fe, 1 = Gd
/// If a species is absent, its outputs are 0.
void compute_bulk_m(const std::vector<uint8_t>& species,
    const std::vector<double>& mx,
    const std::vector<double>& my,
    const std::vector<double>& mz, BulkValues& bulk);

#endif //REDUCTIONS_H
