#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "params.h"

inline void normalize3(double& x, double& y, double& z) {
    if (const double n = std::sqrt(x*x + y*y + z*z); n > 0.0) {
        const double inv = 1.0 / n;
        x *= inv; y *= inv; z *= inv;
    }
}

inline void normalize3(Vec3& vec) {
    normalize3(vec.x, vec.y, vec.z);
}

inline double interpolate_clamped(const std::vector<double>& X,
    const std::vector<double>& Y, const double x)
{
    if (Y.empty())
        throw std::runtime_error("Y should not be empty");
    if (X.size() != Y.size())
        throw std::runtime_error("X and Y should be the same size");
    if (X.size() <= 1) return Y[0];
    if (x <= X.front()) return Y.front();
    if (x >= X.back())  return Y.back();
    const auto it =
        std::upper_bound(X.begin(), X.end(), x);
    const size_t i1 = static_cast<std::size_t>(it - X.begin());
    const size_t i0 = i1 - 1;
    const double w = (x - X[i0]) / (X[i1] - X[i0]);
    return Y[i0] + w * (Y[i1] - Y[i0]);
}

#endif //MATH_UTILS_H
