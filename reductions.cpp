#include "reductions.h"
#include <cmath>

void compute_bulk_m(const std::vector<uint8_t>& species,
    const std::vector<double>& mx,
    const std::vector<double>& my,
    const std::vector<double>& mz, BulkValues& bulk)
{
    double sum_mx_Fe  = 0.0, sum_my_Fe  = 0.0, sum_mz_Fe  = 0.0;
    double sum_mx_Gd  = 0.0, sum_my_Gd  = 0.0, sum_mz_Gd  = 0.0;
    double sum_mx_all = 0.0, sum_my_all = 0.0, sum_mz_all = 0.0;
    int    cnt_Fe = 0, cnt_Gd = 0;

    const int N = static_cast<int>(species.size());
    for (int i = 0; i < N; ++i) {
        const double mx_i = mx[i];
        const double my_i = my[i];
        const double mz_i = mz[i];

        sum_mx_all += mx_i;
        sum_my_all += my_i;
        sum_mz_all += mz_i;

        if (species[i] == 0) // Fe
        {
            sum_mx_Fe += mx_i;
            sum_my_Fe += my_i;
            sum_mz_Fe += mz_i;
            ++cnt_Fe;
        }
        else if (species[i] == 1) // Gd
        {
            sum_mx_Gd += mx_i;
            sum_my_Gd += my_i;
            sum_mz_Gd += mz_i;
            ++cnt_Gd;
        }
    }

    if (N > 0) {
        const double inv = 1.0 / static_cast<double>(N);
        bulk.mx_bulk = sum_mx_all * inv;
        bulk.my_bulk = sum_my_all * inv;
        bulk.mz_bulk = sum_mz_all * inv;
    }
    else {
        bulk.mx_bulk = bulk.my_bulk = bulk.mz_bulk = 0.0;
    }

    if (cnt_Fe > 0) {
        const double inv = 1.0 / static_cast<double>(cnt_Fe);
        bulk.mx_Fe = sum_mx_Fe * inv;
        bulk.my_Fe = sum_my_Fe * inv;
        bulk.mz_Fe = sum_mz_Fe * inv;
    }
    else {
        bulk.mx_Fe = bulk.my_Fe = bulk.mz_Fe = 0.0;
    }

    if (cnt_Gd > 0) {
        const double inv = 1.0 / static_cast<double>(cnt_Gd);
        bulk.mx_Gd = sum_mx_Gd * inv;
        bulk.my_Gd = sum_my_Gd * inv;
        bulk.mz_Gd = sum_mz_Gd * inv;
    }
    else {
        bulk.mx_Gd = bulk.my_Gd = bulk.mz_Gd = 0.0;
    }
}

void compute_bulk_fields(const std::vector<uint8_t>& species,
    const std::vector<double>& Hx_exch_tesla, const std::vector<double>& Hy_exch_tesla, const std::vector<double>& Hz_exch_tesla,
    const std::vector<double>& Hx_anis_tesla, const std::vector<double>& Hy_anis_tesla, const std::vector<double>& Hz_anis_tesla,
    const std::vector<double>& Hx_ther_tesla, const std::vector<double>& Hy_ther_tesla, const std::vector<double>& Hz_ther_tesla,
    BulkFields& bulk_fields)
{
    double sum_Hx_exch_tesla_Fe  = 0.0, sum_Hy_exch_tesla_Fe  = 0.0, sum_Hz_exch_tesla_Fe  = 0.0;
    double sum_Hx_anis_tesla_Fe  = 0.0, sum_Hy_anis_tesla_Fe  = 0.0, sum_Hz_anis_tesla_Fe  = 0.0;
    double sum_Hx_ther_tesla_Fe  = 0.0, sum_Hy_ther_tesla_Fe  = 0.0, sum_Hz_ther_tesla_Fe  = 0.0;

    double sum_Hx_exch_tesla_Gd  = 0.0, sum_Hy_exch_tesla_Gd  = 0.0, sum_Hz_exch_tesla_Gd  = 0.0;
    double sum_Hx_anis_tesla_Gd  = 0.0, sum_Hy_anis_tesla_Gd  = 0.0, sum_Hz_anis_tesla_Gd  = 0.0;
    double sum_Hx_ther_tesla_Gd  = 0.0, sum_Hy_ther_tesla_Gd  = 0.0, sum_Hz_ther_tesla_Gd  = 0.0;

    int cnt_Fe = 0, cnt_Gd = 0;

    const int N = static_cast<int>(species.size());
    for (int i = 0; i < N; ++i) { // 0=Fe, 1=Gd
        if (species[i] == 0) {
            sum_Hx_exch_tesla_Fe += fabs(Hx_exch_tesla[i]);
            sum_Hy_exch_tesla_Fe += fabs(Hy_exch_tesla[i]);
            sum_Hz_exch_tesla_Fe += fabs(Hz_exch_tesla[i]);

            sum_Hx_anis_tesla_Fe += (Hx_anis_tesla[i]);
            sum_Hy_anis_tesla_Fe += (Hy_anis_tesla[i]);
            sum_Hz_anis_tesla_Fe += (Hz_anis_tesla[i]);

            sum_Hx_ther_tesla_Fe += fabs(Hx_ther_tesla[i]);
            sum_Hy_ther_tesla_Fe += fabs(Hy_ther_tesla[i]);
            sum_Hz_ther_tesla_Fe += fabs(Hz_ther_tesla[i]);

            ++cnt_Fe;
        }
        else if (species[i] == 1) {
            sum_Hx_exch_tesla_Gd += fabs(Hx_exch_tesla[i]);
            sum_Hy_exch_tesla_Gd += fabs(Hy_exch_tesla[i]);
            sum_Hz_exch_tesla_Gd += fabs(Hz_exch_tesla[i]);

            sum_Hx_anis_tesla_Gd += (Hx_anis_tesla[i]);
            sum_Hy_anis_tesla_Gd += (Hy_anis_tesla[i]);
            sum_Hz_anis_tesla_Gd += (Hz_anis_tesla[i]);

            sum_Hx_ther_tesla_Gd += fabs(Hx_ther_tesla[i]);
            sum_Hy_ther_tesla_Gd += fabs(Hy_ther_tesla[i]);
            sum_Hz_ther_tesla_Gd += fabs(Hz_ther_tesla[i]);

            ++cnt_Gd;
        }
    }

    if (cnt_Fe > 0) {
        const double inv = 1.0 / static_cast<double>(cnt_Fe);
        bulk_fields.Hx_exch_tesla_Fe = sum_Hx_exch_tesla_Fe * inv;
        bulk_fields.Hy_exch_tesla_Fe = sum_Hy_exch_tesla_Fe * inv;
        bulk_fields.Hz_exch_tesla_Fe = sum_Hz_exch_tesla_Fe * inv;

        bulk_fields.Hx_anis_tesla_Fe = sum_Hx_anis_tesla_Fe * inv;
        bulk_fields.Hy_anis_tesla_Fe = sum_Hy_anis_tesla_Fe * inv;
        bulk_fields.Hz_anis_tesla_Fe = sum_Hz_anis_tesla_Fe * inv;

        bulk_fields.Hx_ther_tesla_Fe = sum_Hx_ther_tesla_Fe * inv;
        bulk_fields.Hy_ther_tesla_Fe = sum_Hy_ther_tesla_Fe * inv;
        bulk_fields.Hz_ther_tesla_Fe = sum_Hz_ther_tesla_Fe * inv;
    }
    else {
        bulk_fields.Hx_exch_tesla_Fe = bulk_fields.Hy_exch_tesla_Fe = bulk_fields.Hz_exch_tesla_Fe = 0.0;
        bulk_fields.Hx_anis_tesla_Fe = bulk_fields.Hy_anis_tesla_Fe = bulk_fields.Hz_anis_tesla_Fe = 0.0;
        bulk_fields.Hx_ther_tesla_Fe = bulk_fields.Hy_ther_tesla_Fe = bulk_fields.Hz_ther_tesla_Fe = 0.0;
    }

    if (cnt_Gd > 0) {
        const double inv = 1.0 / static_cast<double>(cnt_Gd);
        bulk_fields.Hx_exch_tesla_Gd = sum_Hx_exch_tesla_Gd * inv;
        bulk_fields.Hy_exch_tesla_Gd = sum_Hy_exch_tesla_Gd * inv;
        bulk_fields.Hz_exch_tesla_Gd = sum_Hz_exch_tesla_Gd * inv;

        bulk_fields.Hx_anis_tesla_Gd = sum_Hx_anis_tesla_Gd * inv;
        bulk_fields.Hy_anis_tesla_Gd = sum_Hy_anis_tesla_Gd * inv;
        bulk_fields.Hz_anis_tesla_Gd = sum_Hz_anis_tesla_Gd * inv;

        bulk_fields.Hx_ther_tesla_Gd = sum_Hx_ther_tesla_Gd * inv;
        bulk_fields.Hy_ther_tesla_Gd = sum_Hy_ther_tesla_Gd * inv;
        bulk_fields.Hz_ther_tesla_Gd = sum_Hz_ther_tesla_Gd * inv;
    }
    else {
        bulk_fields.Hx_exch_tesla_Gd = bulk_fields.Hy_exch_tesla_Gd = bulk_fields.Hz_exch_tesla_Gd = 0.0;
        bulk_fields.Hx_anis_tesla_Gd = bulk_fields.Hy_anis_tesla_Gd = bulk_fields.Hz_anis_tesla_Gd = 0.0;
        bulk_fields.Hx_ther_tesla_Gd = bulk_fields.Hy_ther_tesla_Gd = bulk_fields.Hz_ther_tesla_Gd = 0.0;
    }
}
