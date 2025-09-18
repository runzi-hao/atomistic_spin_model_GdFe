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

void compute_bulk_fields(const std::vector<uint8_t>& species,
    const std::vector<double> &Hx_exch_tesla, const std::vector<double> &Hy_exch_tesla, const std::vector<double> &Hz_exch_tesla,
    const std::vector<double>& Hx_anis_tesla, const std::vector<double>& Hy_anis_tesla,
    const std::vector<double>& Hz_anis_tesla,
    const std::vector<double>& Hx_ther_tesla, const std::vector<double>& Hy_ther_tesla,
    const std::vector<double>& Hz_ther_tesla,
    BulkFields& bulk_fields);

#endif //REDUCTIONS_H
